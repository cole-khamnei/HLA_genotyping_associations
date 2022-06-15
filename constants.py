import os
import platform
import pathlib

###############################################################################################################
### System Constants ###
###############################################################################################################

if platform.system() == 'Darwin':
    DATA_PATH = "/Users/cole/Documents/_research/rabadan_lab/data"
    dir_divider = "/"
elif platform.system() == 'Windows':
    DATA_PATH = "C://Users//Cole//Documents//_drive//columbia_education//research//rabadan_lab//data"
    dir_divider = "//"

else:
    raise NotImplementedError("Unsupported OS")


###############################################################################################################
### Path Constants ###
###############################################################################################################

HLA_GENOTYPING_DIR_PATH = pathlib.Path(__file__).parent.resolve()

RESOURCES_DIR_PATH = os.path.join(HLA_GENOTYPING_DIR_PATH, "resources")
COVER_PLOTS_DIR_PATH = os.path.join(HLA_GENOTYPING_DIR_PATH, "cover_plots")


UK_BIOBANK_DATA_PATH = DATA_PATH + f"{dir_divider}uk_biobank"

UK_BIOBANK_DEV_DATA_CSV = os.path.join(UK_BIOBANK_DATA_PATH, "small_ukbiobank.csv")
UK_BIOBANK_DATA_CSV = os.path.join(UK_BIOBANK_DATA_PATH, "ukbiobank.csv")

UK_BIOBANK_INDEX_HTML_PATH = os.path.join(UK_BIOBANK_DATA_PATH, "ukbiobank.html")
UK_BIOBANK_INDEX_CSV_PATH = os.path.join(UK_BIOBANK_DATA_PATH, "ukbiobank_index.csv")

UDI_LOOKUP_OUTLIERS_CSV = os.path.join(RESOURCES_DIR_PATH, "udi_lookup_outliers.csv")
UDI_LOOKUP_CORE_CSV = os.path.join(RESOURCES_DIR_PATH, "udi_lookup_core.csv")

###############################################################################################################
### UK BioBank Feature Constants ###
###############################################################################################################

TEMPORARILY_EXCLUDE_FEATURE_TOKENS = ["spirometry_blow_data_points"]

###############################################################################################################
### End ###
###############################################################################################################
