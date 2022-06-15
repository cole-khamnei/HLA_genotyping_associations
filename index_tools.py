
import os
from typing import Any, Dict, Callable, Optional, Tuple

import pandas as pd
import numpy as np
from bs4 import BeautifulSoup as bsoup
from thefuzz import fuzz

import constants

####################################################################################################
### Constants ###
####################################################################################################

UDI_MAPS_ARE_LOADED = False
UDI_TO_NAME_MAP = {}
NAME_TO_UDI_MAP = {}

####################################################################################################
### UK BioBank Index Helpers ###
####################################################################################################


def split_data_coding(data_coding: str) -> Tuple[Optional[int], Optional[int], Optional[str], Optional[str]]:
    """ Splits the UKBB index coding str into parts"""
    if len(data_coding) == 0:
        return [None, None, None, None]

    _, data_coding, _, n_members, data_type, _, _, _, pretype, type_ = data_coding.lower().split()
    return [int(data_coding), int(n_members), data_type, pretype + "_" + type_.rstrip(".")]


def load_index() -> pd.DataFrame:
    """ Loads the UK BioBank index csv and if does not exist, creates it."""
    if os.path.exists(constants.UK_BIOBANK_INDEX_CSV_PATH):
        ukbb_index = pd.read_csv(constants.UK_BIOBANK_INDEX_CSV_PATH)
    else:
        ukbb_html = bsoup(open(constants.UK_BIOBANK_INDEX_HTML_PATH,'r').read())
        ukbb_index_html = ukbb_html.find_all("table")[1]
        ukbb_index = pd.read_html(str(ukbb_index_html))[0]
        ukbb_index.columns = [col.lower() for col in ukbb_index.columns]
        ukbb_index.to_csv(constants.UK_BIOBANK_INDEX_CSV_PATH, index=False)

    ukbb_index["data_coding"] = ukbb_index["description"].apply(lambda desc: desc.split("Uses")[1] if "Uses" in desc else "")
    ukbb_index["description"] = ukbb_index["description"].apply(lambda desc: desc.split("Uses")[0])

    data_coding_info = np.array(ukbb_index["data_coding"].apply(split_data_coding).to_list())
    data_info = pd.DataFrame(data_coding_info, columns=["data_code", "data_n_members", "data_type", "data_structure"])
    ukbb_index = pd.concat([ukbb_index, data_info], axis=1).drop(["data_coding"], axis=1)

    return ukbb_index


def add_biobank_info_to_index(ukbb_index: pd.DataFrame, ukbb_data: pd.DataFrame) -> pd.DataFrame:
    """ Adds relevant statistics from the ukbb data to the index"""

    ukbb_index["counts"] = np.array(ukbb_data.count().tolist())
    ukbb_index["frequency"] = ukbb_index["counts"] / len(ukbb_data)
    return ukbb_index

####################################################################################################
### UK BioBank UDI Helpers ###
####################################################################################################


def load_partial_udi_lookup_map() -> dict:
    """ Loads the UDI lookup tables"""
    core_udi_lookup = pd.read_csv(constants.UDI_LOOKUP_CORE_CSV)
    outlier_udi_lookup = pd.read_csv(constants.UDI_LOOKUP_OUTLIERS_CSV)
    
    partial_udi_lookup = pd.concat([core_udi_lookup, outlier_udi_lookup])

    partial_labeled_udis = partial_udi_lookup.loc[partial_udi_lookup["name"] != "_"]
    partial_udi_to_name_map = dict(zip(partial_labeled_udis["udi"], partial_labeled_udis["name"]))
    return partial_udi_to_name_map


def add_udi_names_to_index(ukbb_index: pd.DataFrame)-> pd.DataFrame:
    """ Adds the names to the index based on udi"""

    partial_udi_to_name_map = load_partial_udi_lookup_map()

    names = []
    for udi in ukbb_index["udi"]:
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

    ukbb_index["name"] = names
    return ukbb_index


class UDIMap:

    def __init__(self, ukbb_index: pd.DataFrame):
        assert "name" in ukbb_index, "'name' column not found in ukbb_index, run add_udi_names_to_index first."

        self.udi_to_name_map = dict(zip(ukbb_index["udi"], ukbb_index["name"]))
        self.name_to_udi_map = dict(zip(ukbb_index["name"], ukbb_index["udi"]))


    def get_udi(self, name: str) -> str:
        """ gets a UDI from a feature name"""

        if isinstance(name, str):
            return self.name_to_udi_map.get(name, name)

        return [self.get_udi(name_i) for name_i in name]


    def get_name(self, udi: str) -> str:
        """ gets the feature name from a udi"""
        
        if isinstance(udi, str):
            return self.udi_to_name_map.get(udi, udi)

        return [self.get_name(udi_i) for udi_i in udi]


    def udi_wrapper(self, function: Callable, *args: Tuple[Any], **kwargs: Dict[str, Any]) -> Any:
        """ wraps a callable function, converting all feature name strings to udi strings."""

        args = {self.get_udi(arg) for arg in args}
        kwargs = {key: self.get_udi(value) if isinstance(value, str) else value for key, value in kwargs.items()}
        return function(*args, **kwargs)


####################################################################################################
### UK BioBank Feature Helpers ###
####################################################################################################


def relevant_feature_search(ukbb_index: pd.DataFrame, term: str) -> pd.DataFrame:
    """ finds  features relevant to the search term."""
    modified_names = ukbb_index["name"].apply(lambda s: s.replace("_", " ") + " " if s else "")
    found_indices = [i for (i, description) in enumerate(modified_names + ukbb_index["description"])
                     if fuzz.partial_ratio(description.lower(), term.lower()) > 95]
    return ukbb_index.iloc[found_indices]


####################################################################################################
### End ###
####################################################################################################
