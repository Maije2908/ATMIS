# Passive element type
#this has to be in a class namespace in order to use it in match/case szenarios
class El:
    INDUCTOR: int = 1
    CAPACITOR: int = 2
    CMC: int = 3

PROMINENCE_DEFAULT = .5 #dB
ACOUSTIC_RESONANCE_PROMINENCE = 1.0 #dB

# offset factor for the first resonance detected after the main peak (min_f = f0 * OFFSET_FACTOR)
# if the first resonance after the main peak is not fitted, consider setting this to a lower value
MIN_ZONE_OFFSET_FACTOR = 2

#multiplication factor for the parasitic element of the main resonance (C for inductor, L for capacitor)
MAIN_RES_PARASITIC_LOWER_BOUND = 0.5
MAIN_RES_PARASITIC_UPPER_BOUND = 2
#max/min values for the main resonance
MIN_R_FE = 10
MAX_R_FE = 1e9
MIN_R_ISO = 1e5
MAX_R_ISO = 1e9
R_ISO_VALUE = 10e6

SERIES_THROUGH = 1
SHUNT_THROUGH = 2



MAX_W_FACTOR = 1.0001
MIN_W_FACTOR = 1/MAX_W_FACTOR
BW_MAX_FACTOR = 1.01
BW_MIN_FACTOR = 1/BW_MAX_FACTOR

# factor to stretch the bandwidth of the last frequency zone (1 = no stretch)
BANDWIDTH_STRETCH_LAST_ZONE = 1


#threshold for the calculation of the offset; necessary for small coils that have a lot of zero crossings at low frequencies
PHASE_OFFSET_THRESHOLD = 60 #60 #°
PHASE_OFFSET_THRESHOLD_CAPS = 20
#value for detection of the inductive/capacitive range; if phase is below this value, inductive/capacitive range will not be detected
PERMITTED_MIN_PHASE = 75 #75


#parameters for the smoothing filter
SAVGOL_WIN_LENGTH_REL = 0.013 # relative savgol window length; yields 52 samples for 4001 point of measurement data
SAVGOL_POLY_ORDER = 2 #polynomial order default:2


#multiplication factor for statistical evaluation of the nominal values; this value will be multiplied to the .50 quanti
#le of the slope and gives the max deviation of the .50 quantile
QUANTILE_MULTIPLICATION_FACTOR = 5

MINIMUM_PRECISION = 1e-12 #if we encounter values that get singular, here is the threshold

DEFAULT_OFFSET_PEAK = 40 #samples; this specifies the default offset for a resonance peak if the 3dB point can't be found

# The relative offset of the data that will be passed to the bandwidth model
BW_MODEL_DATA_OFFSET_STRETCH = 1.5 # (-)

############################## RANGE VARIABLES FOR HIGHER ORDER RESONANCES #############################################

# These variables define how much the higher order resonance parameters are allowed to deviate from the initial guess
# The variables are *factors*, i.e. multiplicative constants. This means a value of RMAX = 2 allows the maximum value of
# a resistor to be two times its initial guess

# resistor ranges for CAPACITORS
RMINFACTOR_CAP = .2
RMAXFACTOR_CAP = 5

# resistor ranges for INDUCTORS
RMINFACTOR_COIL = .5
RMAXFACTOR_COIL = 2

# capacitor ranges for CAPACITORS
CMINFACTOR_CAP = 1e-1
CMAXFACTOR_CAP = 1e1

# capacitor ranges for INDUCTORS
CMINFACTOR_COIL = 1e-1
CMAXFACTOR_COIL = 1e1

############################## ENDS RANGE VARIABLES FOR HIGHER ORDER RESONANCES ########################################



#mode flags
class fcnmode:
    FIT:        int = 1
    FIT_LOG:    int = 6
    OUTPUT:     int = 2
    ANGLE:      int = 3
    FIT_REAL:   int = 4
    FIT_IMAG:   int = 5

class multiple_fit:
    FULL_FIT = 1
    MAIN_RES_FIT = 2

class calc_method:
    SERIES = 1
    SHUNT = 2

class captype:
    GENERIC = 1
    MLCC = 2
    HIGH_C = 3

class cmctype:
    MULTIRESONANCE = 1
    PLATEAU = 2
    NANOCRYSTALLINE = 3

class capunits:
    FARADS:         float = 1
    MILLIFARADS:    float = 1e-3
    MICROFARADS:    float = 1e-6
    NANOFARADS:     float = 1e-9

class indunits:
    HENRYS:         float = 1
    MILLIHENRIES:   float = 1e-3
    MICROHENRIES:   float = 1e-6

class funits:
    HERTZ:          float = 1
    KILOHERTZ:      float = 1e3
    MEGAHERTZ:      float = 1e6

#determines whether to generate differnce plots or not
OUTPUT_DIFFPLOTS = 1

#Debug Plots
DEBUG_BW_MODEL = 0
DEBUG_BW_MODEL_VERBOSE = 0
DEBUG_FIT = 0
DEBUG_MESSAGES = 1
DEBUG_BW_DETECTION = 0
DEBUG_MULTIPLE_FITE_FIT = 0


#determines whether to show bode plots or only save them
SHOW_BODE_PLOTS = False


#logging
LOGGING_VERBOSE = 0

