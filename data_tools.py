import os
import pandas as pd

from typing import Any, Optional, Tuple

import constants, index_tools, medical_code_tools, utilities

if utilities.is_jupyter_notebook():
    from tqdm.notebook import tqdm
else:
    from tqdm import tqdm


########################################################################################################################
### Constants ###
########################################################################################################################


########################################################################################################################
### UK BioBank Data Tools ###
########################################################################################################################


def get_reduced_feature_set_indices(biobank_index: pd.DataFrame) -> pd.Series:
    """ """

    with open(constants.REDUCED_FEATURE_SET_PATH, 'r') as reduced_feature_set_file:
        reduced_feature_set = reduced_feature_set_file.read().split("\n")

    return biobank_index.loc[biobank_index["name"].isin(reduced_feature_set)]["column"].values


def load_biobank_data(dev_mode: bool, udi_map: index_tools.UDIMap, signifier: str = "",
                      biobank_index: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    """ loads the UK BioBank data and converts the udis to common names."""

    data_csv_path = constants.get_uk_biobank_data_csv_path(dev_mode, signifier=signifier)

    if not dev_mode:
        assert biobank_index is not None, "Valid biobank_index must be passed to load_biobank_data if not in dev mode."
        reduced_feature_set_indices = get_reduced_feature_set_indices(biobank_index)

        print("Reduced feature set has", len(reduced_feature_set_indices), "features.")
        biobank_data = pd.read_csv(data_csv_path, low_memory=True, dtype=str, encoding='iso-8859-1',
                                   usecols=reduced_feature_set_indices)

    else:
        biobank_data = pd.read_csv(data_csv_path, low_memory=True, dtype=str, encoding='iso-8859-1')

    biobank_data.columns = udi_map.get_name(biobank_data.columns)
    print(f"UK BioBank Data Loaded.\nSize: {biobank_data.shape[0]} rows x {biobank_data.shape[1]} columns")

    return biobank_data


def clean_biobank_data(biobank_data: pd.DataFrame, biobank_index: pd.DataFrame) -> pd.DataFrame:
    """ cleans the UK BioBank dataframe."""

    biobank_data = biobank_data.loc[:, biobank_data.columns.notna()]
    exclude_features = [column for column in biobank_data.columns
                        if any(token in column for token in constants.TEMPORARILY_EXCLUDE_FEATURE_TOKENS)]
    biobank_data_length = len(biobank_data)
    exclude_features += [feature for feature in biobank_data.columns
                         if biobank_data[feature].isna().sum() == biobank_data_length]

    biobank_data = biobank_data[[feature for feature in biobank_data.columns if feature not in exclude_features]]

    features = biobank_index.query("type == 'Continuous' or type == 'Integer'")["name"]
    return biobank_data.astype({feature: 'float' for feature in features if feature in biobank_data})


def create_reduced_feature_set(biobank_data: pd.DataFrame) -> pd.DataFrame:
    """ builds the a reduced feature set for easier data exploration."""

    reduced_feature_set = {}
    for feature in biobank_data.columns:
        feature_tokens = feature.split("_")
        feature_name_root = "_".join(feature_tokens[:-1]) if "." in feature_tokens[-1] else feature

        if feature < reduced_feature_set.get(feature_name_root, "zzz"):
            reduced_feature_set[feature_name_root] = feature

    return list(reduced_feature_set.values())


def decode_medical_codes(med_code_mapping: medical_code_tools.MedicalCodeMapping,
                         biobank_data: pd.DataFrame) -> pd.DataFrame:
    """ Maps all of the medical codes to values."""
    values = {}
    for feature in tqdm(biobank_data.columns, desc="Mapping Medical Codes", unit=" feature"):
        values[feature] = med_code_mapping.decode(biobank_data[feature], name=feature)

    return pd.DataFrame(values)


def load_all_biobank_components(dev_mode: bool, signifier: str = "") -> Tuple[pd.DataFrame, pd.DataFrame,
                                                                              medical_code_tools.MedicalCodeMapping]:
    """ Returns the biobank data, index, and medical code mapping"""

    print("Importing BioBank Index and Data:")
    with utilities.Timer() as t:
        biobank_index = index_tools.load_index()
        biobank_index = index_tools.add_udi_names_to_index(biobank_index)
        udi_map = index_tools.UDIMap(biobank_index)

        biobank_data = load_biobank_data(dev_mode, udi_map, signifier=signifier, biobank_index=biobank_index)

        biobank_index = index_tools.add_biobank_info_to_index(biobank_index, biobank_data)
        biobank_data = clean_biobank_data(biobank_data, biobank_index)

    if not medical_code_tools.all_biobank_codes_downloaded(biobank_index):
        medical_code_tools.download_all_biobank_codes(biobank_index, overwrite=False, suppress=True)

    med_code_mapping = medical_code_tools.MedicalCodeMapping(biobank_index)

    biobank_data = decode_medical_codes(med_code_mapping, biobank_data)

    return biobank_data, biobank_index, med_code_mapping


def biobank_search(med_code_mapping: medical_code_tools.MedicalCodeMapping,
                   biobank_data: pd.DataFrame, term: str) -> pd.DataFrame:
    """ Searchs for codes that resemble the search term and the counts the occurences in biobank"""
    df = med_code_mapping.search_codes(term).copy(deep=True)
    counts = []

    for feature, meaning in zip(df["name"], df["meaning"]):
        if feature and feature in biobank_data:
            counts.append((biobank_data[feature] == meaning).sum())
        else:
            counts.append(0)
    df["count"] = counts
    return df.sort_values(["count", "code_format"], ascending=False)


########################################################################################################################
### HLA tools ###
########################################################################################################################


def load_hla_data() -> pd.DataFrame:
    """ Loads the HLA allele tsv data"""
    biobank_hla_allele_path = os.path.join(constants.UK_BIOBANK_DATA_PATH,
                                       f"uk_biobank_{constants.UK_BIOBANK_VERSION}_HLA_alleles.tsv")
    hla_alleles = pd.read_csv(biobank_hla_allele_path, sep='\t')
    hla_alleles = hla_alleles[["eid", "A1", "A2", "B1", "B2", "C1", "C2"]].astype({"eid": str})
    
    hla_alleles["zygosity"] = (hla_alleles["A1"] == hla_alleles["A2"]) * 1
    hla_alleles["zygosity"] += (hla_alleles["B1"] == hla_alleles["B2"]) * 1
    hla_alleles["zygosity"] += (hla_alleles["C1"] == hla_alleles["C2"]) * 1
    
    return hla_alleles


########################################################################################################################
### End ###
########################################################################################################################
