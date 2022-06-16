import pandas as pd

import constants, index_tools

from typing import Any, Optional

########################################################################################################################
### Constants ###
########################################################################################################################


########################################################################################################################
### UK BioBank Data Tools ###
########################################################################################################################


def load_biobank_data(data_csv_path: str, udi_map: index_tools.UDIMap) -> pd.DataFrame:
    """ loads the UK BioBank data and converts the udis to common names."""

    biobank_data = pd.read_csv(data_csv_path, low_memory=False)
    biobank_data.columns = udi_map.get_name(biobank_data.columns)
    print(f"UK BioBank Data Loaded.\nSize: {biobank_data.shape[0]} rows x {biobank_data.shape[1]} columns")

    return biobank_data


def clean_biobank_data(biobank_data: pd.DataFrame) -> pd.DataFrame:
    """ cleans the UK BioBank dataframe."""

    biobank_data = biobank_data.loc[:, biobank_data.columns.notna()]
    exclude_features = [column for column in biobank_data.columns
                        if any(token in column for token in constants.TEMPORARILY_EXCLUDE_FEATURE_TOKENS)]
    biobank_data_length = len(biobank_data)
    exclude_features += [feature for feature in biobank_data.columns
                         if biobank_data[feature].isna().sum() == biobank_data_length]

    return biobank_data[[feature for feature in biobank_data.columns if feature not in exclude_features]]


def create_reduced_feature_set(biobank_data: pd.DataFrame) -> pd.DataFrame:
    """ builds the a reduced feature set for easier data exploration."""

    reduced_feature_set = {}
    for feature in biobank_data.columns:
        feature_tokens = feature.split("_")
        feature_name_root = "_".join(feature_tokens[:-1]) if "." in feature_tokens[-1] else feature

        if feature < reduced_feature_set.get(feature_name_root, "zzz"):
            reduced_feature_set[feature_name_root] = feature

    return list(reduced_feature_set.values())


########################################################################################################################
### End ###
########################################################################################################################
