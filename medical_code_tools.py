import os
import requests
import time

import pandas as pd


import utilities

if utilities.is_jupyter_notebook():
    from tqdm.notebook import tqdm
else:
    from tqdm import tqdm

import constants


########################################################################################################################
### Constants ###
########################################################################################################################


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


def download_all_biobank_codes(biobank_index: pd.DataFrame, overwrite: bool = False) -> None:
    """  Downloads all relevant biobank code lookup tables and descriptions."""

    biobank_codes = sorted(biobank_index["data_code"].dropna().unique())
    for code in tqdm(biobank_codes, desc="Downloading BioBank code info", unit=" code"):
        download_biobank_code_data(code, overwrite=overwrite)


########################################################################################################################
### End ###
########################################################################################################################
