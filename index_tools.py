
import os
import requests

from typing import Any, Dict, Callable, Iterable, List, Optional, TypeVar, Tuple, Union

import pandas as pd
import numpy as np

from bs4 import BeautifulSoup as bsoup
from thefuzz import fuzz
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


def split_data_coding(data_coding: str) -> Tuple[Optional[int], Optional[int], Optional[str], Optional[str]]:
    """ Splits the UKBB index coding str into parts"""
    if len(data_coding) == 0:
        return [None, None, None, None]

    _, data_coding, _, n_members, data_type, _, _, _, pretype, type_ = data_coding.lower().split()
    return [int(data_coding), int(n_members), data_type, pretype + "_" + type_.rstrip(".")]


def load_index() -> pd.DataFrame:
    """ Loads the UK BioBank index csv and if does not exist, creates it."""
    if os.path.exists(constants.UK_BIOBANK_INDEX_CSV_PATH):
        biobank_index = pd.read_csv(constants.UK_BIOBANK_INDEX_CSV_PATH)
    else:
        ukbb_html = bsoup(open(constants.UK_BIOBANK_INDEX_HTML_PATH,'r').read())
        biobank_index_html = ukbb_html.find_all("table")[1]
        biobank_index = pd.read_html(str(biobank_index_html))[0]
        biobank_index.columns = [col.lower() for col in biobank_index.columns]
        biobank_index.to_csv(constants.UK_BIOBANK_INDEX_CSV_PATH, index=False)

    biobank_index["data_coding"] = biobank_index["description"].apply(lambda desc: desc.split("Uses")[1] if "Uses" in desc else "")
    biobank_index["description"] = biobank_index["description"].apply(lambda desc: desc.split("Uses")[0])

    data_coding_info = np.array(biobank_index["data_coding"].apply(split_data_coding).to_list())
    data_info = pd.DataFrame(data_coding_info, columns=["data_code", "data_n_members", "data_type", "data_structure"])
    biobank_index = pd.concat([biobank_index, data_info], axis=1).drop(["data_coding"], axis=1)

    return biobank_index


def add_biobank_info_to_index(biobank_index: pd.DataFrame, ukbb_data: pd.DataFrame) -> pd.DataFrame:
    """ Adds relevant statistics from the ukbb data to the index"""

    biobank_index["counts"] = np.array(ukbb_data.count().tolist())
    biobank_index["frequency"] = biobank_index["counts"] / len(ukbb_data)
    return biobank_index

########################################################################################################################
### UK BioBank UDI Helpers ###
########################################################################################################################


def load_partial_udi_lookup_map() -> dict:
    """ Loads the UDI lookup tables"""
    core_udi_lookup = pd.read_csv(constants.UDI_LOOKUP_CORE_CSV)
    outlier_udi_lookup = pd.read_csv(constants.UDI_LOOKUP_OUTLIERS_CSV)
    
    partial_udi_lookup = pd.concat([core_udi_lookup, outlier_udi_lookup])

    partial_labeled_udis = partial_udi_lookup.loc[partial_udi_lookup["name"] != "_"]
    partial_udi_to_name_map = dict(zip(partial_labeled_udis["udi"], partial_labeled_udis["name"]))
    return partial_udi_to_name_map


def add_udi_names_to_index(biobank_index: pd.DataFrame)-> pd.DataFrame:
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
    r_code = requests.post(constants.UK_BIOBANK_CODING_DOWNLOAD_URL, data={"id":3}).text
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
    #     download_biobank_coding_data(code, overwrite=False)
        pass


########################################################################################################################
### End ###
########################################################################################################################
