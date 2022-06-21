import pandas as pd

import constants, index_tools

from typing import Any, Optional, Tuple

import medical_code_tools, utilities

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


def load_biobank_data(data_csv_path: str, udi_map: index_tools.UDIMap) -> pd.DataFrame:
    """ loads the UK BioBank data and converts the udis to common names."""

    biobank_data = pd.read_csv(data_csv_path, low_memory=False, dtype=str)
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


def decode_medical_codes(med_code_mapping: medical_code_tools.MedicalCodeMapping, biobank_data: pd.DataFrame) -> pd.DataFrame:
    """ Maps all of the medical codes to values."""
    values = {}
    for feature in tqdm(biobank_data.columns, desc="Mapping ICD10 Codes", unit=" feature"):
        values[feature] = med_code_mapping.decode(biobank_data[feature], name=feature)
    
    return pd.DataFrame(values)


def load_all_biobank_components() -> Tuple[pd.DataFrame, pd.DataFrame, medical_code_tools.MedicalCodeMapping]:
    """ Returns the biobank data, index, and medical code mapping"""
  
    print("Importing BioBank Index and Data:")
    with utilities.Timer() as t:
        biobank_index = index_tools.load_index()
        biobank_index = index_tools.add_udi_names_to_index(biobank_index)
        udi_map = index_tools.UDIMap(biobank_index)

        biobank_data = load_biobank_data(constants.UK_BIOBANK_DATA_CSV_PATH, udi_map)

        biobank_index = index_tools.add_biobank_info_to_index(biobank_index, biobank_data)
        biobank_data = clean_biobank_data(biobank_data, biobank_index)

    if not medical_code_tools.all_biobank_codes_downloaded(biobank_index):
        medical_code_tools.download_all_biobank_codes(biobank_index, overwrite=False, suppress=True)
        
    med_code_mapping = medical_code_tools.MedicalCodeMapping(biobank_index)

    biobank_data = decode_medical_codes(med_code_mapping, biobank_data)
    
    return biobank_data, biobank_index, med_code_mapping


def biobank_search(med_code_mapping: medical_code_tools.MedicalCodeMapping, biobank_data: pd.DataFrame, term: str) -> pd.DataFrame:
    """ Searchs for codes that resemble the search term and the counts the occurences in biobank"""
    df = med_code_mapping.search_codes(term).copy(deep=True)
    counts = []
    
    for feature, meaning in zip(df["name"], df["meaning"]):
        if feature:
            counts.append((biobank_data[feature] == meaning).sum())
        else:
            counts.append(0)
    df["count"] = counts
    return df.sort_values(["count", "code_format"], ascending=False)


########################################################################################################################
### End ###
########################################################################################################################
