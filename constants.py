import os
import platform
import pathlib

########################################################################################################################
### System Constants ###
########################################################################################################################

if platform.system() == 'Darwin':
    DATA_PATH = "/Users/cole/Documents/_research/rabadan_lab/data"
    dir_divider = "/"
elif platform.system() == 'Windows':
    DATA_PATH = "C://Users//Cole//Documents//_drive//columbia_education//research//rabadan_lab//data"
    dir_divider = "//"

else:
    raise NotImplementedError("Unsupported OS")

########################################################################################################################
### Path Constants ###
########################################################################################################################

HLA_GENOTYPING_DIR_PATH = pathlib.Path(__file__).parent.resolve()

RESOURCES_DIR_PATH = os.path.join(HLA_GENOTYPING_DIR_PATH, "resources")
COVER_PLOTS_DIR_PATH = os.path.join(HLA_GENOTYPING_DIR_PATH, "cover_plots")
COVER_PLOTS_GENERIC_FILE_PATH = os.path.join(COVER_PLOTS_DIR_PATH, "{}")

UK_BIOBANK_DATA_PATH = DATA_PATH + f"{dir_divider}uk_biobank"

UK_BIOBANK_DEV_DATA_CSV = os.path.join(UK_BIOBANK_DATA_PATH, "small_ukbiobank.csv")
UK_BIOBANK_DATA_CSV = os.path.join(UK_BIOBANK_DATA_PATH, "ukbiobank.csv")

UK_BIOBANK_INDEX_HTML_PATH = os.path.join(UK_BIOBANK_DATA_PATH, "ukbiobank.html")
UK_BIOBANK_INDEX_CSV_PATH = os.path.join(UK_BIOBANK_DATA_PATH, "ukbiobank_index.csv")

UDI_LOOKUP_OUTLIERS_CSV = os.path.join(RESOURCES_DIR_PATH, "udi_lookup_outliers.csv")
UDI_LOOKUP_CORE_CSV = os.path.join(RESOURCES_DIR_PATH, "udi_lookup_core.csv")

CODING_INFO_DIR_PATH = os.path.join(RESOURCES_DIR_PATH, "code_info")

########################################################################################################################
### URL Path Constants ###
########################################################################################################################

UK_BIOBANK_CODING_URL = "https://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id={code}"
UK_BIOBANK_CODING_DOWNLOAD_URL = "https://biobank.ndph.ox.ac.uk/showcase/codown.cgi"

NDPH_DATABASE_ERROR_TOKEN = "<!--\nDATABASE ERROR"

########################################################################################################################
### UK BioBank Feature Constants ###
########################################################################################################################

TEMPORARILY_EXCLUDE_FEATURE_TOKENS = ["spirometry_blow_data_points"]

########################################################################################################################
### End ###
########################################################################################################################
