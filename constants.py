import os
import platform
import pathlib

########################################################################################################################
### System Constants ###
########################################################################################################################

if platform.system() == 'Darwin':
    DATA_PATH = "/Users/cole/Documents/_research/rabadan_lab/data"
    GRANTHAM_DISTANCE_PATH = "/Users/cole/Documents/_research/rabadan_lab/"
    dir_divider = "/"
elif platform.system() == 'Windows':
    DATA_PATH = "C://Users//Cole//Documents//_drive//columbia_education//research//rabadan_lab//data"
    dir_divider = "//"
    GRANTHAM_DISTANCE_PATH = None

else:
    raise NotImplementedError("Unsupported OS")

########################################################################################################################
### Version Constants (EDIT these) ###
########################################################################################################################

# UK_BIOBANK_VERSION = "24512"
UK_BIOBANK_VERSION = "31127"
DEV_MODE = True

UK_BIOBANK_FINALIZED_VERSIONS = ["24512", "31127"]

########################################################################################################################
### Path Constants ###
########################################################################################################################

HLA_GENOTYPING_DIR_PATH = pathlib.Path(__file__).parent.resolve()

RESOURCES_DIR_PATH = os.path.join(HLA_GENOTYPING_DIR_PATH, "resources")
COVER_PLOTS_DIR_PATH = os.path.join(HLA_GENOTYPING_DIR_PATH, "cover_plots")
COVER_PLOTS_GENERIC_FILE_PATH = os.path.join(COVER_PLOTS_DIR_PATH, "{}")

UK_BIOBANK_DATA_PATH = DATA_PATH + f"{dir_divider}uk_biobank_{UK_BIOBANK_VERSION}"

UK_BIOBANK_DEV_DATA_CSV = os.path.join(UK_BIOBANK_DATA_PATH, f"uk_biobank_{UK_BIOBANK_VERSION}_small.csv")
UK_BIOBANK_DATA_CSV = os.path.join(UK_BIOBANK_DATA_PATH, f"uk_biobank_{UK_BIOBANK_VERSION}.csv")

UK_BIOBANK_INDEX_HTML_PATH = os.path.join(UK_BIOBANK_DATA_PATH, f"uk_biobank_{UK_BIOBANK_VERSION}.html")
UK_BIOBANK_INDEX_CSV_PATH = os.path.join(UK_BIOBANK_DATA_PATH, f"uk_biobank_{UK_BIOBANK_VERSION}_index.csv")

UK_BIOBANK_HLA_ALLELES_TSV_PATH = os.path.join(UK_BIOBANK_DATA_PATH, f"uk_biobank_{UK_BIOBANK_VERSION}_HLA_alleles.tsv")
UK_BIOBANK_HLA_ALLELES_TSV_FULL_PATH = os.path.join(UK_BIOBANK_DATA_PATH, f"uk_biobank_{UK_BIOBANK_VERSION}_HLA_alleles_full.tsv")

UDI_LOOKUP_OUTLIERS_CSV = os.path.join(RESOURCES_DIR_PATH, "udi_lookup_outliers.csv")
UDI_LOOKUP_CORE_CSV = os.path.join(RESOURCES_DIR_PATH, "udi_lookup_core.csv")

UDI_LOOKUP_DIR_PATH = os.path.join(RESOURCES_DIR_PATH, "udi_lookups")
UDI_LOOKUP_VERSION_GENERIC_PATH = os.path.join(UDI_LOOKUP_DIR_PATH, "udi_lookup_{version}.csv")

REDUCED_FEATURE_SET_PATH = os.path.join(RESOURCES_DIR_PATH, "feature_resources", "reduced_feature_set.txt")

CODING_INFO_DIR_PATH = os.path.join(RESOURCES_DIR_PATH, "code_info")

########################################################################################################################
### Mode Set Constants ###
########################################################################################################################

# UK_BIOBANK_DATA_CSV_PATH = UK_BIOBANK_DEV_DATA_CSV if DEV_MODE else UK_BIOBANK_DATA_CSV


def get_uk_biobank_data_csv_path(dev_mode: bool, signifier: str = "") -> str:
    """ Returns the UK_BIOBANK_DATA_CSV_PATH based on the current dev mode."""
    path = UK_BIOBANK_DEV_DATA_CSV if dev_mode else UK_BIOBANK_DATA_CSV
    return path.replace(".csv", f"_{signifier}.csv") if signifier != "" else path


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
### Variable Constants ###
########################################################################################################################

DEFAULT = "default"

########################################################################################################################
### End ###
########################################################################################################################
