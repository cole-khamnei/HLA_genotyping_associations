import os
import glob
import requests
import time

from collections.abc import Iterable
from typing import Optional, TypeVar, Union

import pandas as pd
import treelib as tl

import utilities

if utilities.is_jupyter_notebook():
    from tqdm.notebook import tqdm
else:
    from tqdm import tqdm

import constants



########################################################################################################################
### Constants ###
########################################################################################################################

S = TypeVar('S')
ArrayOrItem = Union[Iterable[S], S]

########################################################################################################################
### UK BioBank Coding Lookup formatting.###
########################################################################################################################


def download_biobank_code_data(code: int, overwrite: bool = False, suppress: bool = False) -> None:
    """ Downloads the coding information for a data code."""

    lookup_path = os.path.join(constants.CODING_INFO_DIR_PATH, f"code_lookup_{code}.csv")
    description_path = os.path.join(constants.CODING_INFO_DIR_PATH, f"code_description_{code}.txt")

    if os.path.exists(lookup_path) and os.path.exists(description_path) and not overwrite:
        if not suppress:
            print(f"code {code} files already exist.")
        return

    # Parse the description to provide variable description info
    r_desc = requests.get(constants.UK_BIOBANK_CODING_URL.format(code=code)).text
    assert not r_desc.startswith(constants.NDPH_DATABASE_ERROR_TOKEN), "NDPH database is down. Please wait."

    # Parse description here.
    code_description = r_desc.split("<br>Name: ")[1]
    name, code_description = code_description.split("\n<p>Description: ")
    description, data_format = code_description.split("\n<p>")[:2]

    # Download coding information and format as csv.
    r_code = requests.post(constants.UK_BIOBANK_CODING_DOWNLOAD_URL, data={"id": code}).text
    time.sleep(1)
    assert not r_code.startswith(constants.NDPH_DATABASE_ERROR_TOKEN), "NDPH database is down. Please wait."
    r_code = r_code.replace('\t', ',')

    with open(lookup_path, 'w') as code_lookup_file:
        code_lookup_file.write(r_code)

    with open(description_path, 'w') as code_description_file:
    	code_description_file.write(name + "\n")
    	code_description_file.write(description + "\n")
    	code_description_file.write(data_format + "\n")


def download_all_biobank_codes(biobank_index: pd.DataFrame, overwrite: bool = False, suppress: bool = False) -> None:
    """  Downloads all relevant biobank code lookup tables and descriptions."""

    biobank_codes = sorted(biobank_index["data_code"].dropna().unique())
    for code in tqdm(biobank_codes, desc="Downloading BioBank code info", unit=" code"):
        download_biobank_code_data(code, overwrite=overwrite, suppress=suppress)


########################################################################################################################
### Medical Code Mapping###
########################################################################################################################


class MedicalCodeMapping:
    def __init__(self, biobank_index: pd.DataFrame):
        self.name_to_code_format_map = dict(biobank_index[["name", "data_code"]].values)
        self.code_format_to_name_map = dict(biobank_index[["data_code", "name"]].values)
        self.code_format_lookup = self.get_code_format_info_lookup()

        self.code_formats_with_tree_structure = []
        self.ICD10_master_list = []
        for code_format, info in self.code_format_lookup.items():
            if "node_id" in info["coded_values_df"]:
                self.code_formats_with_tree_structure.append(code_format)
            ICD10_code_meanings = info["coded_values_df"][["coding","meaning"]].copy(deep=True)
            ICD10_code_meanings["code_format"] = code_format
            names = biobank_index.query(f"data_code == '{code_format}'").sort_values("name")["name"]
            ICD10_code_meanings["feature"] = names.values[0] if len(names) > 0 else None
            self.ICD10_master_list.append(ICD10_code_meanings)

        self.ICD10_master_list = pd.concat(self.ICD10_master_list).drop_duplicates(["coding", "meaning", "code_format"])

    def clean_name(self, name: str) -> str:
        """ Cleans the name"""
        return name
        if "." not in name:
            return name

        return "_".join(name.split("_")[:-1])
        
    def get_code_format_from_name(self, name: str) -> str:
        """ Gets the medical code format used by variable 'name'"""
        name = self.clean_name(name)
        return self.name_to_code_format_map[name]
    
    def get_code_lookup_from_name(self, name: str):
        """ Gets the medical code format lookup from variable 'name'"""
        name = self.clean_name(name)
        return self.code_format_lookup[self.get_code_format_from_name(name)]
    
    def name_has_code_format(self, name: str) -> bool:
        """ Checks if a variable 'name' has a code format"""
        name = self.clean_name(name)
        return name in self.name_to_code_format_map

    def parse_code_format(self, csv_path: str) -> pd.DataFrame:
        """ """
        with open(csv_path, 'r') as csv_file:
            csv_text = csv_file.read()

        lines = csv_text.split("\n")
        header, lines = lines[0].split(","), lines[1:]
        meaning_index = header.index("meaning")
        assert meaning_index == 1, "Meaning index in code lookup csv is not the second column."

        values = []
        for line in lines:
            tokens = line.split(",")
            if len(header) > 2:
                end_index = 2 - len(header) 
                values.append(tokens[:1] + [";".join(tokens[1:end_index])] + tokens[end_index:])
            else:
                values.append(tokens[:1] + [";".join(tokens[1:])])

        return pd.DataFrame(values[1:], columns=header)

    def build_tree(self, name: Optional[str] = None, code_format: Optional[str] = None) -> tl.Tree:
        """ """
        assert bool(code_format) != bool(name), "Only name or code should be provided."
        if name:
            name = self.clean_name(name)
            code = self.get_code_format_from_name(name)
        else:
            name = self.code_format_to_name_map.get(code_format, code_format)
            code = code_format


        tree = tl.Tree()
        tree.create_node(name, "0")

        df = self.code_format_lookup[code]["coded_values_df"]
        items = [row for i, row in df.sort_values("node_id").iterrows()]
        index = 0

        while len(items) > 0:
            row = items[index]
            if tree.contains(row["parent_id"]):
                tree.create_node(row["meaning"] + f" ({row['coding']})", row["node_id"], parent=row["parent_id"])
                del items[index]
                index = 0
                fail_stop = False
            else:
                if index < len(items) - 1:
                    index = index + 1
                elif fail_stop:
                    break
                else:
                    fail_stop = True
                    index = 0

        tree.show()
        return tree

    def get_code_format_info_lookup(self) -> dict:
        """ Creates a code format info lookup dict."""

        code_format_info_lookup = {}
        for code_format_csv_path in glob.glob(os.path.join(constants.CODING_INFO_DIR_PATH, "code_lookup_*.csv")):
            code_format = code_format_csv_path.split("code_lookup_")[1].rstrip(".csv")
            code_format_description_path = code_format_csv_path.replace("lookup",
                                                                        "description").replace(".csv", ".txt")

            with open(code_format_description_path, 'r') as code_format_description_file:
                name, description, data_format = code_format_description_file.read().split("\n", maxsplit=2)

            code_values_df = self.parse_code_format(code_format_csv_path)
            code_values = dict(zip(code_values_df["coding"], code_values_df["meaning"]))

            code_format_info_lookup[code_format] = {"name": name, "desc": description,
                                             "data_format": data_format, "coded_values": code_values,
                                             "coded_values_df": code_values_df}
        return code_format_info_lookup
    
    def decode(self, coded_value: str, code_format: Optional[str] = None, name: Optional[str] = None):
        """ Decodes a coded value"""
        assert bool(code_format) != bool(name), "Only name or code should be provided."
        if name:
            name = self.clean_name(name)
        
        code_format = self.get_code_format_from_name(name) if name else code_format
        if code_format is None:
            return coded_value
        
        if isinstance(coded_value, Iterable) and not isinstance(coded_value, str):
            return [self.decode(coded_value_i, code_format=code_format) for coded_value_i in coded_value]
        
        return self.code_format_lookup[code_format]["coded_values"].get(str(coded_value), coded_value)


    def search_codes(self, search_terms: ArrayOrItem[str], fuzzy_threshold: int = 95, and_search: bool = False) -> pd.DataFrame:
        """ Searches the descriptions of ICD codes for relevent features"""
        indices = utilities.fuzzy_index_search(search_terms, self.ICD10_master_list["meaning"].values, fuzzy_threshold=fuzzy_threshold, and_search=and_search)
        return self.ICD10_master_list.iloc[indices]


########################################################################################################################
### End ###
########################################################################################################################
