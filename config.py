
import constants

FULL_FIT = True
MULTIPROCESSING_COUNT = 16

FREQ_UPPER_LIMIT = 2e9
FREQ_LOWER_LIMIT = 0

MAX_ORDER = 15

FORCE_SINGLE_POINT_MODEL = False

CAPUNIT = constants.capunits.NANOFARADS
INDUNIT = constants.indunits.MICROHENRIES
FUNIT   = constants.funits.MEGAHERTZ

R_FILL_CAP = 10e6
R_FILL_IND = 1e-1

FIT_BY = constants.fcnmode.FIT


CMC_REQUIRED_CONFIGURATIONS = ["DM", "CM"]

#threshold for the calculation of the offset; necessary for small coils that have a lot of zero crossings at low frequencies
PHASE_OFFSET_THRESHOLD = 60 #°
#value for detection of the inductive/capacitive range; if phase is below this value, inductive/capacitive range will not be detected
PERMITTED_MIN_PHASE = 75 #°
