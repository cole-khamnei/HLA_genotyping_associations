
import os
from typing import Any, Dioct, Callable, Optional, Tuple

import pandas as pd
import numpy as np
from bs4 import BeautifulSoup as bsoup

import constants

####################################################################################################
### Constants ###
####################################################################################################



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


####################################################################################################
### UK BioBank UDI Helpers ###
####################################################################################################


def get_udi(name: str) -> str:
    """ """
    if isinstance(name, str):
        return name_to_udi_map.get(name, name)

    return [get_udi(name_i) for name_i in name]


def get_name_from_udi(udi: str) -> str:
    """ """
    if isinstance(udi, str):
        return udi_to_name_map.get(udi, udi)

    return [get_name_from_udi(udi_i) for udi_i in udi]


def udi_wrapper(function: Callable, *args: Tuple[Any], **kwargs: Dict[str, Any]) -> Any:
    """ """
    args = {get_udi(arg) for arg in args}
    kwargs = {key: get_udi(value) if isinstance(value, str) else value for key, value in kwargs.items()}
    return function(*args, **kwargs)


def relevant_feature_search(ukbb_index: pd.DataFrame, term: str) -> pd.DataFrame:
    """ """
    modified_names = ukbb_index["name"].apply(lambda s: s.replace("_", " ") + " " if s else "")
    found_indices = [i for (i, description) in enumerate(modified_names + ukbb_index["description"])
                     if fuzz.partial_ratio(description.lower(), term.lower()) > 95]
    return ukbb_index.iloc[found_indices]


####################################################################################################
### End ###
####################################################################################################
