
import argparse
import os
import requests
import glob
import time

from typing import Any, Dict, Callable, Iterable, List, Optional, TypeVar, Tuple, Union

import pandas as pd
import numpy as np

from bs4 import BeautifulSoup as bsoup
from thefuzz import fuzz

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
### UK BioBank Index Helpers ###
########################################################################################################################


def overwrite_override(file: Optional[str] = ""):
    """ Checks the version of biobank used and prevents overwrites of any files that are *finalized*"""
    if constants.UK_BIOBANK_VERSION in constants.UK_BIOBANK_FINALIZED_VERSIONS:
        assert False, f"Overwrite attempt of file {file} for version {constants.UK_BIOBANK_VERSION} which is finalized."


########################################################################################################################
### UK BioBank Index Helpers ###
########################################################################################################################


def create_biobank_index(biobank_html_path: str, overwrite: bool = False) -> None:
    """ Creates a csv with the index information for the biobank."""

    biobank_index_path = biobank_html_path.replace(".html", "_index.csv")
    if os.path.exists(biobank_index_path) and not overwrite:
        print(f"\n{biobank_index_path} already exists.\nIt can be overwritten with the 'overwrite' argument")
        return

    biobank_html = bsoup(open(biobank_html_path, 'r').read(), features="lxml")
    biobank_index_html = biobank_html.find_all("table")[1]
    biobank_index = pd.read_html(str(biobank_index_html))[0]
    biobank_index.columns = [col.lower() for col in biobank_index.columns]
    biobank_index.to_csv(biobank_index_path, index=False)
    print(f"\nSuccessflly created {biobank_index_path}")


def extract_data_coding(description: str) -> str:
    """ Extracts the data code from the description."""
    if "data-coding" not in description:
        return None

    return description.split("data-coding ")[1].split()[0]


def load_index() -> pd.DataFrame:
    """ Loads the UK BioBank index csv and if does not exist, creates it."""

    biobank_index = pd.read_csv(constants.UK_BIOBANK_INDEX_CSV_PATH)

    biobank_index["data_code"] = biobank_index["description"].apply(extract_data_coding)
    biobank_index["description"] = biobank_index["description"].apply(lambda desc: desc.split("Uses")[0])

    return biobank_index


def add_biobank_info_to_index(biobank_index: pd.DataFrame, ukbb_data: pd.DataFrame) -> pd.DataFrame:
    """ Adds relevant statistics from the ukbb data to the index"""

    biobank_index["counts"] = np.array(ukbb_data.count().tolist())
    biobank_index["frequency"] = biobank_index["counts"] / len(ukbb_data)
    return biobank_index

########################################################################################################################
### UK BioBank UDI Lookup Creators ###
########################################################################################################################


def get_indices_missing_names(biobank_index: pd.DataFrame) -> pd.DataFrame:
    """ Finds all the indices that are missing a name variable."""

    missing_name_index = biobank_index.loc[biobank_index["name"].isna()].copy(deep=True)
    missing_name_index["primary_udi"] = missing_name_index["udi"].apply(lambda s: s.split("-")[0])
    missing_name_index = missing_name_index.drop_duplicates(["primary_udi"])
    missing_name_index["name"] = "_"
    udi_lookup_columns = ["column", "count", "type", "description", "udi", "name"]
    return missing_name_index[udi_lookup_columns]


def create_udi_lookup_file(version: str, missing_name_index: pd.DataFrame, overwrite: bool = False) -> None:
    """ Writes the udi lookup file to be filled out."""

    if overwrite:
        overwrite_override()

    udi_lookup_file = constants.UDI_LOOKUP_VERSION_GENERIC_PATH.format(version=version)

    if os.path.exists(udi_lookup_file):
        print(f"Version {version} UDI lookup file already exists.")
        if not overwrite:
            print("If you want to overwrite, use 'overwrite=True'")
            return

        print(f"Overwriting existing version {version} UDI lookup file.")
    missing_name_index.to_csv(udi_lookup_file, index=False)
    print(f"Successfully created UDI lookup file:\n{udi_lookup_file}")


########################################################################################################################
### UK BioBank UDI Helpers ###
########################################################################################################################


def load_partial_udi_lookup_map() -> dict:
    """ Loads the UDI lookup tables"""

    udi_lookup_paths = glob.glob(os.path.join(constants.UDI_LOOKUP_DIR_PATH, "udi_lookup_*.csv"))
    partial_udi_lookup = pd.concat([pd.read_csv(udi_lookup_path) for udi_lookup_path in udi_lookup_paths])

    partial_labeled_udis = partial_udi_lookup.loc[partial_udi_lookup["name"] != "_"]
    partial_udi_to_name_map = dict(zip(partial_labeled_udis["udi"], partial_labeled_udis["name"]))
    return partial_udi_to_name_map


def add_udi_names_to_index(biobank_index: pd.DataFrame) -> pd.DataFrame:
    """ Adds the names to the index based on udi"""

    partial_udi_to_name_map = load_partial_udi_lookup_map()

    names = []
    for udi in biobank_index["udi"]:
        if "-" not in udi or udi.endswith("-0.0"):
            names.append(partial_udi_to_name_map.get(udi, None))
        else:
            udi_stem, udi_modifier = udi.split("-")
            new_name = None

            if udi_stem + "-0.0" in partial_udi_to_name_map:
                names.append(f"{partial_udi_to_name_map[udi_stem + '-0.0']}_{udi_modifier}")
            elif udi_stem + "-0.1" in partial_udi_to_name_map:
                names.append(f"{partial_udi_to_name_map[udi_stem + '-0.1']}_{udi_modifier}")
            else:
                names.append(None)

    biobank_index["name"] = names
    return biobank_index


class UDIMap:

    def __init__(self, biobank_index: pd.DataFrame):
        assert "name" in biobank_index, "'name' column not found in biobank_index, run add_udi_names_to_index first."

        self.udi_to_name_map = dict(zip(biobank_index["udi"], biobank_index["name"]))
        self.name_to_udi_map = dict(zip(biobank_index["name"], biobank_index["udi"]))

    def get_udi(self, name: ArrayOrItem[str]) -> ArrayOrItem[str]:
        """ gets a UDI from a feature name"""

        if isinstance(name, str):
            return self.name_to_udi_map.get(name, name)

        return [self.get_udi(name_i) for name_i in name]

    def get_name(self, udi: ArrayOrItem[str]) -> ArrayOrItem[str]:
        """ gets the feature name from a udi"""

        if isinstance(udi, str):
            return self.udi_to_name_map.get(udi, udi)

        return [self.get_name(udi_i) for udi_i in udi]

    def udi_wrapper(self, function: Callable, *args: Tuple[Any], **kwargs: Dict[str, Any]) -> Any:
        """ wraps a callable function, converting all feature name strings to udi strings."""

        args = {self.get_udi(arg) for arg in args}
        kwargs = {key: self.get_udi(value) if isinstance(value, str) else value for key, value in kwargs.items()}
        return function(*args, **kwargs)


########################################################################################################################
### UK BioBank Feature Helpers ###
########################################################################################################################


def fuzzy_index_search(term: str, descriptions: Iterable[str], fuzzy_threshold: int = 95,
                       and_search: bool = False) -> List[int]:
    """ Searches a list of descriptions and returns any with a fuzzy index above a threshold."""

    if isinstance(term, str):
        term = term.lower()
        return [i for (i, desc) in enumerate(descriptions) if fuzz.partial_ratio(desc.lower(), term) > fuzzy_threshold]

    index_sets = [set(fuzzy_index_search(term_i, descriptions, fuzzy_threshold=fuzzy_threshold)) for term_i in term]

    final_indices = index_sets[0]
    for indices in index_sets[1:]:
        if and_search:
            final_indices = final_indices.intersection(indices)
        else:
            final_indices = final_indices.union(indices)

    return sorted(final_indices)


def term_search(biobank_index: pd.DataFrame, search_terms: ArrayOrItem[str],
                fuzzy_threshold: int = 95, and_search: bool = False) -> pd.DataFrame:
    """ searches the biobank_index for relevant features."""

    search_terms = [search_terms] if isinstance(search_terms, str) else search_terms
    descriptions = biobank_index["name"].apply(lambda s: s.replace("_", " ") + " " if s else "")
    descriptions = descriptions + biobank_index["description"]

    indices = fuzzy_index_search(search_terms, descriptions, fuzzy_threshold=fuzzy_threshold, and_search=and_search)
    return biobank_index["name"].iloc[indices].tolist()


########################################################################################################################
### UK BioBank Coding Lookup formatting.###
########################################################################################################################


def download_biobank_code_data(code: int, overwrite: bool = False) -> None:
    """ Downloads the coding information for a data code."""

    lookup_path = os.path.join(constants.CODING_INFO_DIR_PATH, f"code_lookup_{code}.csv")
    description_path = os.path.join(constants.CODING_INFO_DIR_PATH, f"code_description_{code}.txt")

    if os.path.exists(lookup_path) and os.path.exists(description_path) and not overwrite:
        print(f"code {code} files already exist.")
        return

    # Parse the description to provide variable description info
    r_desc = requests.get(constants.UK_BIOBANK_CODING_URL.format(code=code)).text
    assert not r_desc.startswith(constants.NDPH_DATABASE_ERROR_TOKEN), "NDPH database is down. Please wait."

    # Parse description here.
    code_description = r_desc

    # Download coding information and format as csv.
    r_code = requests.post(constants.UK_BIOBANK_CODING_DOWNLOAD_URL, data={"id": code}).text
    time.sleep(1)
    assert not r_code.startswith(constants.NDPH_DATABASE_ERROR_TOKEN), "NDPH database is down. Please wait."
    r_code = r_code.replace('\t', ',')

    with open(lookup_path, 'w') as code_lookup_file:
        code_lookup_file.write(r_code)

    with open(description_path, 'w') as code_description_file:
        code_description_file.write(r_desc)


def download_all_biobank_codes(biobank_index: pd.DataFrame, overwrite: bool = False) -> None:
    """  Downloads all relevant biobank code lookup tables and descriptions."""

    biobank_codes = sorted(biobank_index["data_code"].dropna().unique())
    for code in tqdm(biobank_codes, desc="Downloading BioBank code info", unit=" code"):
        download_biobank_code_data(code, overwrite=overwrite)


########################################################################################################################
### Main ###
########################################################################################################################


def argument_parser() -> str:
    """"""
    parser = argparse.ArgumentParser(description='Index biobank csv')
    parser.add_argument('--biobank', dest="biobank_path", help='Path to biobank csv',
                        type=str, required=True)
    parser.add_argument('--overwrite', dest="overwrite", action="store_true",
                        help='overwrites existing index / UDI lookup csv.')
    parser.add_argument('--create-index', dest="create_index", action="store_true", help='create index.')
    parser.add_argument('--create-udi-lookup', dest="create_udi_lookup", action="store_true", help='create UDI lookup.')
    parser.add_argument('--full', dest="full_prep", action="store_true", help='create index and UDI lookup.')

    inputs = parser.parse_args()

    assert inputs.biobank_path.endswith(".csv"), f"Biobank path must be a csv, given:'{inputs.biobank_path}'."
    assert os.path.exists(inputs.biobank_path), f"BioBank CSV '{inputs.biobank_path}' does not exist."

    if inputs.full_prep:
        mode = "full"
    elif inputs.create_udi_lookup:
        mode = "udi_lookup_creation"
    elif inputs.create_index:
        mode = "index_creation"
    else:
        raise ValueError("No mode specified use flags '--full', '--create-index', or '--create-udi-lookup'")

    return inputs.biobank_path, mode, inputs.overwrite


def main():
    biobank_path, mode, overwrite = argument_parser()

    if mode in ["full", "index_creation"]:
        biobank_html_path = biobank_path.replace(".csv", ".html")
        create_biobank_index(biobank_html_path, overwrite=overwrite)

    if mode in ["full", "udi_lookup_creation"]:
        biobank_index_full = load_index()
        biobank_index_full = add_udi_names_to_index(biobank_index_full)

        missing_name_index = get_indices_missing_names(biobank_index_full)
        if len(missing_name_index) > 0:
            create_udi_lookup_file(constants.UK_BIOBANK_VERSION, missing_name_index, overwrite=overwrite)


if __name__ == '__main__':
    main()

########################################################################################################################
### End ###
########################################################################################################################
