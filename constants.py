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
    GRANTHAM_DISTANCE_PATH = "C://Users//Cole//Documents//_drive//columbia_education//research//rabadan_lab//"

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
### Lists ###
########################################################################################################################

FEMALE_SPECIFIC_CANCERS = ["breast cancer", "female genital tract cancer", "cervical cancer",
                           "cin/pre-cancer cells cervix", "fallopian tube cancer", "ovarian cancer",
                           "uterine/endometrial cancer", "vaginal cancer", "vulval cancer" ]

MALE_SPECIFIC_CANCERS = ["male genital tract cancer", "penis cancer", "prostate cancer", "testicular cancer"]
SEX_SPECIFIC_CANCERS = FEMALE_SPECIFIC_CANCERS + MALE_SPECIFIC_CANCERS

ICD10_CANCER_TYPES = {
    "C00":" lip",
    "C01":" tongue base",
    "C02":" tongue other",
    "C03":" gum",
    "C04":" mouth floor",
    "C05":" mouth palate",
    "C06":" mouth other",
    "C07":" parotid gland",
    "C08":" salivary glands other",
    "C09":" tonsil",
    "C10":" oropharynx",
    "C11":" nasopharynx",
    "C12":" pyriform sinus",
    "C13":" hypopharynx",
    "C15":" esophagus",
    "C16":" stomach",
    "C17":" small intestine",
    "C18":" colon",
    "C19":" rectosigmoid junction",
    "C20":" rectal",
    "C21":" anal",
    "C22":" liver",
    "C23":" gallbladder",
    "C24":" biliary tract other",
    "C25":" pancreas",
    "C26":" digestive organs other",
    "C30":" nasal cavity and middle ear",
    "C31":" accessory sinuses",
    "C32":" larynx",
    "C33":" trachea",
    "C34":" lung",
    "C37":" thymus",
    "C38":" heart, mediastinum and pleura",
    "C43":" melanoma",
    "C44":" non melanoma skin cancer",
    "C45":" mesothelioma",
    "C46":" kaposi's sarcoma",
    "C50":" breast",
    "C51":" vulva",
    "C52":" vagina",
    "C53":" cervix uteri",
    "C54":" corpus uteri",
    "C55":" uterus other",
    "C56":" ovary",
    "C57":" female genital other",
    "C58":" placenta",
    "C60":" penis",
    "C61":" prostate",
    "C62":" testis",
    "C63":" male genital other",
    "C64":" kidney",
    "C65":" renal pelvis",
    "C66":" ureter",
    "C67":" bladder",
    "C68":" urinary other",
    "C69":" eye and adnexa",
    "C70":" meninges",
    "C71":" brain",
    "C72":" cns other",
    "C73":" thyroid gland",
    "C74":" adrenal gland",
    "C75":" endocrine other",
    "C76":" other primary",
    "C77":" lymph nodes secondary",
    "C78":" respiratory and digestive organs secondary",
    "C79":" other secondary",
    "C81":" hodgkin's disease",
    "C82":" follicular [nodular] non-hodgkin's lymphoma",
    "C83":" diffuse non-hodgkin's lymphoma",
    "C84":" peripheral and cutaneous t-cell lymphomas",
    "C85":" non-hodgkin's lymphoma other",
    "C86":" t/nk-cell lymphoma other",
    "C88":" malignant immunoproliferative diseases",
    "C90":" multiple myeloma",
    "C91":" lymphoid leukaemia",
    "C92":" myeloid leukaemia",
    "C93":" monocytic leukaemia",
    "C94":" leukaemia other",
    "C95":" leukaemia unspecified",
    "D00":" carcinoma in situ of oral cavity, oesophagus and stomach155",
    "D02":" carcinoma in situ of middle ear and respiratory system101",
    "D03":" melanoma in situ",
    "D04":" carcinoma in situ of skin",
    "D05":" carcinoma in situ of breast",
    "D06":" carcinoma in situ of cervix uteri"
}

########################################################################################################################
### End ###
########################################################################################################################
