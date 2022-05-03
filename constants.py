import platform

###############################################################################################################
### Constants ###
###############################################################################################################

if platform.system() == 'Darwin':
    DATA_PATH = "/Users/cole/Documents/_research/rabadan_lab/data"
    dir_divider = "/"
elif platform.system() == 'Windows':
    DATA_PATH = "C://Users//Cole//Documents//_drive//columbia_education//research//rabadan_lab//data"
    dir_divider = "//"

else:
    raise NotImplementedError("Unsupported OS")        

UK_BIOBANK_DATA_PATH = DATA_PATH + f"{dir_divider}uk_biobank"

###############################################################################################################
### End ###
###############################################################################################################
