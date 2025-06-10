import lmfit
import matplotlib.pyplot as plt
import numpy
import numpy as np
import scipy
import skrf
from scipy import signal
from lmfit import minimize, Parameters
from scipy.signal import find_peaks
import decimal
import copy
import sys
import logging


# import constants
import constants
from config import MAX_ORDER
from constants import *
import config



class Fitter:
    """
    TODO: docstring; this will be used when help() command is invoked
    """

    def __init__(self, name, freq, data, fit_type, shunt_series = SERIES_THROUGH, captype = captype.GENERIC,
                 logger_instance = logging.getLogger(), nominal_value = None, series_resistance = None,
                 peak_detection_prominence = PROMINENCE_DEFAULT):

        if len(freq) != len(data):
            raise Exception("Frequency vector and data vector do not have the same size")

        # Set arguments that need to be passed directly
        self.ser_shunt = shunt_series
        self.fit_type = fit_type
        self.captype = captype
        self.logger = logger_instance

        if peak_detection_prominence is None:
            self.prominence = PROMINENCE_DEFAULT
        else:
            self.prominence = peak_detection_prominence

        self.name = name

        # Set frequency vector and data
        self.freq = freq[(freq > config.FREQ_LOWER_LIMIT) & (freq < config.FREQ_UPPER_LIMIT)]
        self.z21_data = data[(freq > config.FREQ_LOWER_LIMIT) & (freq < config.FREQ_UPPER_LIMIT)]

        # Smooth data
        savgol_length = int(np.floor(SAVGOL_WIN_LENGTH_REL* len(freq))) if int(np.floor(SAVGOL_WIN_LENGTH_REL* len(freq))) > 2 else 3
        self._smooth_data(savgol_length, SAVGOL_POLY_ORDER)

        # Calculate the linear range offset except for High_C model
        if captype != constants.captype.HIGH_C:
            self._calculate_linear_range_offset()
        else:
            self._offset = 0

        # Calculate nominal value if it's not provided
        # if High_C model, nominal value needs to be provided, so calculate_nominal_value should not be invoked
        if nominal_value is None:
            self.calculate_nominal_value()
        else:
            self.get_main_resonance()
            self.nominal_value = nominal_value

        # Calculate series resistance, if not provided
        if series_resistance is None:
            self.calculate_nominal_Rs()
        else:
            self.series_resistance = series_resistance

        # Create an lmfit.Parameters() object to store parameters in
        self.parameters = Parameters()

        # Set fit order to 0 as to have a failsafe if the resonance detection fails
        self.order = 0

        #self.acoustic_resonance_frequency = None

    ########################### CONSTRUCTORS ###########################################################################

    @classmethod
    def from_s2p_file(cls, file, fit_type, shunt_series = SERIES_THROUGH, captype = captype.GENERIC,
                 logger_instance = logging.getLogger(), Z0=50, nominal_value = None, series_resistance = None,
                 peak_detection_prominence = PROMINENCE_DEFAULT):

        # Calculate impedance data
        if shunt_series == SERIES_THROUGH:
            data = cls.calc_series_thru(file, Z0)
        elif shunt_series == SHUNT_THROUGH:
            data = cls.calc_shunt_thru(file, Z0)

        # Get frequency vector
        freq = file.f
        # Get name of file
        name = file.name

        return cls(name, freq, data, fit_type, shunt_series = shunt_series, captype = captype,
                   logger_instance = logger_instance, nominal_value = nominal_value,
                   series_resistance = series_resistance,peak_detection_prominence = peak_detection_prominence)

    ####################################################################################################################
    # Pre-Processing Methods
    ####################################################################################################################

    @staticmethod
    def calc_series_thru(file: skrf.Network, Z0: int = 50) -> numpy.ndarray:
        """
        Function to calculate the series-through impedance of the DUT based on S_21

        :param file: a sNpfile of type skrf.Network()
        :param Z0: (optional) nominal impedance of the measurement system; default is 50
        :return: a vector containing the series through impedance of the DUT (complex ndarray)
        """

        z21_data = 2 * Z0 * ((1 - file.s[:, 1, 0]) / file.s[:, 1, 0])
        return z21_data

    @staticmethod
    def calc_shunt_thru(file: skrf.Network, Z0: int = 50) -> numpy.ndarray:
        """
        Function to calculate the shunt-through impedance of the DUT based on S_21

        :param file: a sNpfile of type skrf.Network()
        :param Z0: (optional) nominal impedance of the measurement system; default is 50
        :return: a vector containing the shunt through impedance of the DUT (complex ndarray)
        """

        z21_data = (Z0 * file.s[:, 1, 0]) / (2 * (1 - file.s[:, 1, 0]))
        return z21_data

    def _smooth_data(self, window, poly_order):
        """
        Function to smooth the impedance data

        A Savitzky-Golay filter is used to smooth the impedance data. The smoothed data will be stored in instance
        variables **data_mag** and **data_ang** for magnitude and phase respectively

        :return: None (stores smoothed data in instance variables data_mag and data_ang)
        """

        sav_gol_mode = 'interp'
        self.data_mag = signal.savgol_filter(abs(self.z21_data), window, poly_order, mode=sav_gol_mode)
        self.data_ang = signal.savgol_filter(np.angle(self.z21_data, deg=True), window,
                                                   poly_order, mode=sav_gol_mode)
        #limit the phase data to +/- 90°
        self.data_ang = np.clip(self.data_ang, -90, 90)

    def calculate_nominal_value(self):
        """
        Function to calculate the nominal value of the DUT, if it was not provided.

        This works by looking for the linear range of the passive element, i.e. the range where the DUT behaves like a
        linear coil/cap. The value is then calculated by using the mean of all obtained values of the linear range

        :return:            Nominal value of the DUT in Henry/Farad (the value is also written to an instance variable)
        :raises Exception:  If the phase of the dataset is too low (i.e. inductive/capacitive linear range can not be
            detected); there is a constant that can be modified in the GUI_config to allow lower phase, however this might make
            the calculation less precise
        """

        try:
            self.get_main_resonance()
        except:
            raise

        #TODO: this part kind of violates the DRY paradigm; don't know how to do that any other way tho
        # this is a late task, we can think about this when everything else works
        match self.fit_type:
            case El.INDUCTOR:

                # Check if the phase of the dataset has a valid range; if the phase is not around 90° we cannot assume
                # inductive range -> throw exception in this case
                if max(self.data_ang[self._offset:self._f0_index]) < PERMITTED_MIN_PHASE:
                    raise Exception("Error: Inductive range not detected (max phase = {value}°).\n"
                                    "Please specify nominal inductance.".format(value=np.round(max(self.data_ang), 1)))

                # Crop data to [offset:f0] in order to obtain the linear range for the calculation of nominal value
                curve_data = self.z21_data[self.freq < self.f0][self._offset:]
                w_data = (self.freq[self.freq < self.f0][self._offset:])*2*np.pi

                # create an array filled with possible values for L; calculation is L = imag(Z)/omega
                L_vals = []
                for it, curve_sample in enumerate(zip(curve_data, w_data)):
                    L_vals.append(np.imag(curve_sample[0])/curve_sample[1])

                # calculate the slope of the magnitude and get the 50% quantile of it; after that find the max slope
                slope_quantile_50 = np.quantile(np.gradient(self.data_mag)[self.freq<self.f0],0.5)
                max_slope = slope_quantile_50 * QUANTILE_MULTIPLICATION_FACTOR
                
                #boolean index the data that has lower than max slope and calculate the mean of it
                L_vals_mean = np.array(L_vals)[np.gradient(self.data_mag)[self.freq<self.f0][self._offset:] < max_slope]
                
                # finally write the obtained nominal value to the instance variable; also write back the offset
                self.nominal_value = np.median(L_vals_mean)

                #TODO: logging etc
                output_dec = decimal.Decimal("{value:.3E}".format(value=self.nominal_value))
                self.logger.info("Nominal Inductance calculated: " + output_dec.to_eng_string())


            case El.CAPACITOR:

                if min(self.data_ang[self._offset:self._f0_index]) > -PERMITTED_MIN_PHASE:
                    raise Exception("Error: Capacitive range not detected (min phase = {value}°).\n"
                                    "Please specify nominal capacitance.".format(value=np.round(min(self.data_ang), 1)))


                #crop data to [offset:f0] in order to find the linear range for the calculation of nominal value
                curve_data = self.z21_data[self.freq < self.f0][self._offset:]
                w_data = (self.freq[self.freq < self.f0][self._offset:])*2*np.pi

                #create an array filled with possible values for C; calculation is C= -1/(imag(Z)*omega)
                C_vals = []
                for it, curve_sample in enumerate(zip(curve_data, w_data)):
                    C_vals.append(-1/(np.imag(curve_sample[0])*curve_sample[1]))

                #TODO: old detection method, check if new one works better/worse/equal

                # # find the 50% quantile of the slope data and define the max slope allowed
                # slope_quantile_50 = np.quantile(np.gradient(C_vals), 0.5)
                # max_slope = slope_quantile_50 * QUANTILE_MULTIPLICATION_FACTOR
                # # boolean index the data that has lower than max slope and calculate the mean
                # C_vals_eff = np.array(C_vals)[abs(np.gradient(C_vals)) < abs(max_slope)]

                #TODO: this is the new detection

                # calculate the slope of the magnitude and get the 50% quantile of it; after that find the max slope
                slope_quantile_50 = np.quantile(np.gradient(self.data_mag)[self.freq < self.f0], 0.5)
                max_slope = slope_quantile_50 * QUANTILE_MULTIPLICATION_FACTOR

                # boolean index the data that has lower than max slope and calculate the mean of it
                C_vals_eff = np.array(C_vals)[np.gradient(self.data_mag)[self.freq < self.f0][self._offset:] < max_slope]

                self.nominal_value = np.median(C_vals_eff)

                output_dec = decimal.Decimal("{value:.3E}".format(value=self.nominal_value))
                self.logger.info("Nominal Capacitance calculated: " + output_dec.to_eng_string())

    def calculate_nominal_Rs(self):
        """
        Function to calculate the series resistance of the DUT if it was not provided.

        Here the minimum of the impedance data is taken; this is **very** imprecise for inductors but works very well
        for capacitors
        :return: obtained series resistance (also writes the value to an instance variable)
        """

        R_s_input = min(abs(self.data_mag))
        self.series_resistance = R_s_input
        output_dec = decimal.Decimal("{value:.3E}".format(value=R_s_input))
        self.logger.info("Nominal Resistance not provided, calculated: " + output_dec.to_eng_string())

        return self.series_resistance

    def _calculate_linear_range_offset(self):
        """
        Method to calculate the offset for the linear range.

        The offset is the number of samples where the phase of the dataset does not yet display inductive/capacitive
        behaviour (i.e. 90°/-90° phase respectively). In order for the main resonant frequency and nominal
        value detection to work, proper inductive/capacitive behaviour is needed.


        :raises IndexError: If the threshold did not get crossed by any value of the phase data
        :raises AttributeError: If the method was invoked before the fit_type instance variable was set.
            Or if the phase data is not present.

        :return: the number of samples to skip in order to be in the linear range safely
            ("safe" in this context means a phase lower/higher than the one specified in #TODO: add constant to description)
        :rtype: int

        """

        # Find the offset for the linear range by taking the first point where the phase is greater/smaller than the
        # specified minimum/maximum for inductors/capacitors respectively
        try:
            match self.fit_type:
                case El.INDUCTOR:
                        offset = np.argwhere(self.data_ang > constants.PHASE_OFFSET_THRESHOLD)[0][0]
                case El.CAPACITOR:
                        offset = np.argwhere(self.data_ang < -constants.PHASE_OFFSET_THRESHOLD)[0][0]

        # If such a crossing of the threshold hasn't been found, raise an exception
        except IndexError:

            # Provide some information about the min/max phase
            match self.fit_type:
                case El.INDUCTOR:
                    minmaxstring = "Max Phase = " + str(np.round(max(self.data_ang),1)) + "°"
                case El.CAPACITOR:
                    minmaxstring = "Min Phase = " + str(np.round(min(self.data_ang),1)) + "°"

            # Raise the exception
            raise IndexError("Error: Bad phase data. "+minmaxstring+"\n"
                            "Determination of the linear inductive/capacitive range not possible.\n"
                            "You can try to lower the "
                            f"phase threshold (currently [+/-]{PHASE_OFFSET_THRESHOLD}°)")

        # If all went well write offset to instance variable
        else:
            self._offset = offset

    def get_main_resonance(self) -> float:
        """
        Method to calculate the main resonant frequency of the DUT.

        This works by looking for the zero crossing of the phase.

        :return: resonant frequency f0 in Hz (also writes the value and the index of the zero crossing to instance variables)
        :rtype: float
        :raises Exception: if the resonant frequency could not be determined
        """

        # Check if offset calculation has been done and if not perform offset calculation
        if self._offset is None:
            self.logger.warn("Offset calculation has not been performed. Performing offset calculation: ")
            # Pass exception if one was caught by the offset detection
            try:
                self._calculate_linear_range_offset()
            except Exception:
                self.logger.info("Failed")
                raise
            else:
                self.logger.info("Success")

        # Set f0 to 0 in order to have feedback, if the method didn't work
        f0 = 0

        # Find zero crossing of the phase
        match self.fit_type:
            case El.INDUCTOR:
                try:
                    index_angle_smaller_zero = np.argwhere(self.data_ang[self._offset:] < 0)
                    index_ang_zero_crossing = self._offset + index_angle_smaller_zero[0][0]
                except Exception:
                    raise Exception("Could not find phase zero crossing for main resonance detection!")
                # TODO: continuity check should not be performed by taking 10 samples after the f0, seems kind of arbitrary
                continuity_check = index_angle_smaller_zero[10][0]

            case El.CAPACITOR:
                try:
                    index_angle_larger_zero = np.argwhere(self.data_ang[self._offset:] > 0)
                    index_ang_zero_crossing = index_angle_larger_zero[0][0] + self._offset
                except Exception:
                    raise Exception("Could not find phase zero crossing for main resonance detection!")
                #TODO: see above
                continuity_check = index_angle_larger_zero[10][0]

        if continuity_check:
            f0 = self.freq[index_ang_zero_crossing]
            self.f0 = f0
            self._f0_index = index_ang_zero_crossing


            output_dec = decimal.Decimal("{value:.3E}".format(value=f0))
            self.logger.info("Detected f0: "+ output_dec.to_eng_string())
            print("Detected f0: "+output_dec.to_eng_string())


        if f0 == 0:
            raise Exception('ERROR: Main resonant frequency could not be determined.')

        return f0

    def get_resonances(self):
        """
        Method to get the higher order resonances.

        Using scipy.find_peaks on the magnitude data, resonant peaks are obtained. Also the bandwidth will be detected
        and a flag will be set if the bandwidth could not be found (due to missing -3dB points)

        :return:    bandwidth list-list of the detected resonances containing lower, center and upper frequencies
            (upper and lower frequency are at the -3dB points); also stores the obtained bands, the peak heights, the order
            (i.e. total number of peaks) and a flag if the -3dB points are not present in instance variables
        :raises Exception: if no resonances have been found
        """

        #TODO: this method should be refactored, but it's a pain

        # create one figure for debug plots
        if constants.DEBUG_BW_DETECTION:
            plt.figure()
            plt.title(self.name)


        freq = self.freq
        magnitude_data = self.data_mag
        phase_data = self.data_ang

        # frequency limit the data
        magnitude_data = magnitude_data[freq < config.FREQ_UPPER_LIMIT]
        freq = freq[freq < config.FREQ_UPPER_LIMIT]

        # get prominence
        prominence_mag = self.prominence

        # find peaks of Magnitude Impedance curve (using scipy.signal.find_peaks)
        # using dB here in order to match the prominence values
        match self.fit_type:
            case constants.El.INDUCTOR:
                mag_maxima = find_peaks(20*np.log10(magnitude_data), prominence=prominence_mag)
            case constants.El.CAPACITOR:
                mag_maxima = find_peaks(-20 * np.log10(magnitude_data), prominence=prominence_mag)


        # get frequency of acquired magnitude maxima
        f_mag_maxima = freq[mag_maxima[0]]

        #ignore all peaks that lie "before" the main resonance and that are to close to the main resonance
        min_zone_start = self.f0 * constants.MIN_ZONE_OFFSET_FACTOR

        mag_maxima_pos = f_mag_maxima[f_mag_maxima > min_zone_start]
        mag_maxima_index = mag_maxima[0][f_mag_maxima > min_zone_start]


        number_zones = len(mag_maxima_pos)
        bandwidth_list = []
        peak_heights = []
        bad_BW_flag = np.zeros((number_zones,2))

        for num_maximum in range(0, number_zones):
            # get index, frequency and value of impedance at that point
            res_fq = mag_maxima_pos[num_maximum]
            res_index = mag_maxima_index[num_maximum]
            res_value = magnitude_data[res_index]

            #get 3dB value
            match self.fit_type:
                case constants.El.INDUCTOR:
                    bw_value = res_value / np.sqrt(2)
                case constants.El.CAPACITOR:
                    bw_value = res_value*np.sqrt(2)

            try:
                #find the index where the 3db value is reached; also check if the frequency is lower than the resonance,
                #but higher than the min zone; if that does not work use the default offset
                #NOTE: since we need the first value in front of the resonance we have to flipud the array
                match self.fit_type:
                    case constants.El.INDUCTOR:
                        f_lower_index = np.flipud(np.argwhere(np.logical_and(freq > min_zone_start, np.logical_and(freq < res_fq, (magnitude_data) < (bw_value)))))[0][0]
                    case constants.El.CAPACITOR:
                        f_lower_index = np.flipud(np.argwhere(np.logical_and(freq > min_zone_start, np.logical_and(freq < res_fq, (magnitude_data) > (bw_value)))))[0][0]
            except IndexError:
                f_lower_index = res_index - constants.DEFAULT_OFFSET_PEAK
                bad_BW_flag[num_maximum][0] = 1

            try:
                match self.fit_type:
                    case constants.El.INDUCTOR:
                        f_upper_index = np.argwhere(np.logical_and(freq > res_fq, (magnitude_data) < (bw_value)))[0][0]
                    case constants.El.CAPACITOR:
                        f_upper_index = np.argwhere(np.logical_and(freq > res_fq, (magnitude_data) > (bw_value)))[0][0]
            except IndexError:
                #here we need to account for the fact that we could overshoot the max index
                if res_index + constants.DEFAULT_OFFSET_PEAK < len(freq):
                    f_upper_index = res_index + constants.DEFAULT_OFFSET_PEAK
                    bad_BW_flag[num_maximum][1] = 1
                else:
                    f_upper_index = len(freq) - 1

            # check if the found 3dB points are in an acceptable range i.e. not "behind" the next peak or "in front of"
            # the previous peak. If that is the case we set the index to a default offset to get a "bandwidth"
            if num_maximum != 0:
                if f_lower_index < mag_maxima_index[num_maximum - 1]:
                    f_lower_index = res_index - constants.DEFAULT_OFFSET_PEAK
                    bad_BW_flag[num_maximum] = 1
            if num_maximum < number_zones-1:
                if f_upper_index > mag_maxima_index[num_maximum + 1]:
                    #again we could overshoot the max index here
                    if res_index + constants.DEFAULT_OFFSET_PEAK < len(freq):
                        f_upper_index = res_index + constants.DEFAULT_OFFSET_PEAK
                        bad_BW_flag[num_maximum][1] = 1
                    else:
                        f_upper_index = len(freq) - 1
                        bad_BW_flag[num_maximum][1] = 1

            f_tuple = [freq[f_lower_index], res_fq, freq[f_upper_index]]
            bandwidth_list.append(f_tuple)
            peak_heights.append(abs(res_value))


        try:
            # here the last resonance can be "streched"; can be useful to "fill" the end zone of the curve
            # NOTE: the stretch factor is set to 1 most of the time, so there is no stretching happening
            stretch_factor = constants.BANDWIDTH_STRETCH_LAST_ZONE
            bandwidth_list[-1][2] = bandwidth_list[-1][2] * stretch_factor
            bandwidth_list[-1][0] = bandwidth_list[-1][0] * stretch_factor

        except IndexError:
            self.logger.info("INFO: No resonances found except the main resonance, consider a lower value for the prominence")

        # THIS IS FOR TESTING
        if constants.DEBUG_BW_DETECTION:
            markerson=np.where(np.isin(self.freq, bandwidth_list))
            plt.loglog(self.freq, self.data_mag, '-bD', markevery=list(markerson[0]))


        self.peak_heights = peak_heights
        self.bandwidths = bandwidth_list
        self.bad_bandwidth_flag = bad_BW_flag
        self.order = len(self.bandwidths)

        return bandwidth_list

    def fit_acoustic_resonance(self, param_set: lmfit.Parameters = None):
        """
        Method to fit the acoustic resonance of MLCCs


        :param param_set: a Parameters() object; should contain the parameters for the main resonance but not for higher order circuits

        :returns: a Parameters() object that is a copy of the Parameters() object passed to the function but containing
            the fitted circuit elements for the acoustic resonance (this Parameters() object is also stored in an instance variable)

        """
        freq = self.freq
        data = self.z21_data

        if param_set is None:
            param_set = self.parameters


        try:
            self.acoustic_resonance_frequency
        except AttributeError as e:
            raise type(e)(str(e) +'. Make sure to set the acoustic resonance frequency '
                                  'before invoking this method!').with_traceback(sys.exc_info()[2])

        # frequency limit the data to frequencies lower than the main resonance
        mag_data_lim = self.data_mag[self.freq < self.f0]
        freq_lim = self.freq[self.freq < self.f0]

        # Find first index of data where we are above the acoustic resonance
        res_fq = self.acoustic_resonance_frequency
        res_index = np.argwhere(self.freq > res_fq)[0][0]

        # Get the value of the data at this point and calculate the approximate value for the -3dB points
        res_value = data[res_index]
        bw_value = abs(res_value) * np.sqrt(2)

        # find upper 3dB point
        try:
            f_upper_index = np.argwhere(np.logical_and(freq > res_fq, (self.data_mag) > (bw_value)))[0][0]
            fu=freq[f_upper_index]
        except IndexError:
            # here we need to account for the fact that we could overshoot the max index
            if res_index + constants.DEFAULT_OFFSET_PEAK < len(freq):
                f_upper_index = res_index + int(constants.DEFAULT_OFFSET_PEAK / 2)
            else:
                f_upper_index = len(freq_lim) - 1

        # find lower 3dB point
        try:
            f_lower_index = np.flipud(np.argwhere((np.logical_and(freq < res_fq, (self.data_mag) > (bw_value)))))[0][0]
            fl = freq[f_lower_index]
        except IndexError:
            f_lower_index = res_index - int(constants.DEFAULT_OFFSET_PEAK / 2)

        freq_mdl = freq_lim[f_lower_index-10:f_upper_index+10]
        data_mdl = mag_data_lim[f_lower_index-10:f_upper_index+10]

        # use the bandwidth model to get estimates for the RLC circuit
        [bl,bu,R,L,C] = self.model_bandwidth(freq_mdl, data_mdl, res_fq)

        # correct the effect the main resonance has on the peak height
        main_res_here = self._calculate_Z(param_set, res_fq, 2, 0, 1, constants.fcnmode.OUTPUT)
        data_here = data[freq<=res_fq][0]
        w_c = res_fq * 2 * np.pi
        Q = res_fq / (bu - bl)

        R_new = abs(1 / (1 / data_here - 1 / main_res_here))
        C_new = 1/(R_new*w_c*Q)

        C=C_new
        params = copy.copy(param_set)

        #scale parameters to match units
        C = C / config.CAPUNIT
        L = L / config.INDUNIT
        w_c = w_c / config.FUNIT
        expr_string_L = '(1/((C_A*'+ str(config.CAPUNIT)+')*(w_A*'+str(config.FUNIT)+')**2))/'+str(config.INDUNIT)

        #create the parameters for the acoustic resonance
        params.add('R_A', value=R_new, min=R*0.8, max=R_new*1.2, vary=True)
        params.add('C_A', value=C, min=C * 0.8, max=C * 1.5, vary=True)
        params.add('w_A', value =w_c, min = w_c*0.9, max=w_c*1.2, vary=True)
        params.add('L_A', value=L, expr=expr_string_L)


        #frequency limit the data to the bandwidth of the circuit and do a fit using the limited data
        modelfreq = freq[np.logical_and(freq > bl, freq < bu)]
        modeldata = data[np.logical_and(freq > bl, freq < bu)]

        out1 = minimize(self._calculate_Z, params,
                        args=(modelfreq, modeldata, 0, 0, config.FIT_BY,),
                        method='powell', options={'xtol': 1e-18, 'disp': True})

        # copy the parameters of the fit result
        params = out1.params

        # finally set the parameters to not vary anymore
        self.change_parameter(params,'R_A',vary=False)
        self.change_parameter(params,'L_A',vary=False)
        self.change_parameter(params,'C_A',vary=False)
        self.change_parameter(params,'w_A',vary=False)

        # write to instance variable and return obtained parameters
        self.parameters = params
        return params

    def get_acoustic_resonance(self, offset = 100000.0) -> float:
        """
        Method to detect an acoustic resonance.
        Finds a peak in the data that has lower frequency than the main resonance.

        :param offset: Frequency offset from the main resonance. This can help avoid misdetection of the main resonance
            as the acoustic resonance. Default is 100 kHz.

        :return: Obtained resonant frequency in Hz; can be None, if no acoustic resonance has been found; also sets the
                instance variable *acoustic_resonance_frequency*
        """

        # Limit the data and frequency to frequencies lower than the main resonant frequency (f0)
        mag_data_lim = self.data_mag[self.freq < self.f0-offset]
        freq_lim = self.freq[self.freq < self.f0-offset]

        # Find any peaks in the data
        mag_maxima = find_peaks(-20 * np.log10(mag_data_lim), height=-200, prominence=ACOUSTIC_RESONANCE_PROMINENCE)

        # Find the peak with the highest prominence
        index_ac_res = np.argwhere(mag_maxima[1]['prominences'] == max(mag_maxima[1]['prominences']))[0][0]

        try:
            res_index = mag_maxima[0][index_ac_res]
            res_fq = freq_lim[res_index]
        except:
            # If we didn't find an acoustic resonance simply set to None
            res_fq = None

        self.acoustic_resonance_frequency = res_fq
        return res_fq

    def create_hi_C_parameters(self, param_set: lmfit.Parameters = None) -> lmfit.Parameters:
        """
        Function to create the parameter set for the bathtub model (electrolytic capacitor model).
        Note: if an additional resonance is present, this method will try to find estimates for the parameters of that
        resonance as well. If it is not pronounced enough, it will only create parameters for the base bathtub model

        :param param_set: (optional) if you want to provide a custom parameter set, parameters will be added to that;
            should be an empty set. If not provided instance variable 'self.parameters' will be used

        :return: parameter set with the parameters for the bathtub model. Also writes parameter set to instance variable
            'self.parameters

        """
        # If no custom parameter set has been supplied, take the instances parameters
        if param_set is None:
            param_set = self.parameters

        C = self.nominal_value / config.CAPUNIT
        param_set.add('C', value = C, min = C*1e-1, max = C*1e1, vary = False)
        param_set.add('R_iso', value = constants.R_ISO_VALUE, min = constants.MIN_R_ISO, max = constants.MAX_R_ISO, vary = False)
        param_set.add('R_s', value = self.series_resistance, min = self.series_resistance*1e-1, max = self.series_resistance*1e1, vary = False)

        # If a resonance has been detected, calculate the value of the inductance of the main resonant circuit to not
        # get lower than 10 times the impedance of the first resonance
        if any(self.bandwidths):
            R1 = self.peak_heights[0]
            f1 = self.bandwidths[0][1]
            w1 = 2 * np.pi * f1
            Z_min = R1 * 1.0

            # Values for the main resonance
            L_min = (Z_min/w1) / config.INDUNIT
            L_value = L_min * 1e1
            L_max = L_min * 1e2
            param_set.add('L', value = L_value, min = L_min, max = L_max)


            # Find estimates for the first resonance

            # Get indices of the band
            b_l = self.bandwidths[0][0]
            b_u = self.bandwidths[0][2]
            f_l_index = np.argwhere(self.freq >= b_l)[0][0]
            f_u_index = np.argwhere(self.freq <= b_u)[-1][0]

            # Get data for bandwidth model
            freq_BW_mdl = self.freq[f_l_index:f_u_index]
            data_BW_mdl = self.data_mag[f_l_index:f_u_index] * np.exp(1j * np.radians(self.data_ang[f_l_index:f_u_index]))

            # Try to estimate the values for the additional resonance via the bandwidth model;
            # sometimes this procedure fails because the resonance is not pronounced enough, in this case we skip it
            # because the basic bathtub model will suffice for modeling
            try:
                [bl, bu, R1, L1, C1] = self.model_bandwidth(freq_BW_mdl,data_BW_mdl, f1)

                L1 = L1 / config.INDUNIT
                C1 = C1 / config.CAPUNIT
                w1 = w1 / config.FUNIT

                param_set.add('R1', value = R1, min = R1*0.5, max = R1*2, vary =True )
                param_set.add('C1', value=C1, min=C1 * 0.3, max=C1 * 3, vary = True)
                param_set.add('w1', value=w1, vary = False)
                param_set.add('BW1', value=((bu-bl))/config.FUNIT, vary = False)
                param_set.add('L1', expr = '(1/( '+ '((w1 * '+str(config.FUNIT)+') ** 2)' + '*(C1*' + str(config.CAPUNIT) + ')) )/'+str(config.INDUNIT)  )

                self.order = 1

            # If the bandwidth model failed, it is very likely that the detected resonance was in fact not a resonance
            # but rather a misdetection, need to determine inductance via dataset
            except:
                # Reset order to 0, so we don't get problems later and calculate the value of the inductance
                self.order = 0
                L = self.calc_L_electrolytic_cap(self.freq, self.z21_data)

                # Write to parameter set with scaling
                L_scaled = L / config.INDUNIT
                param_set.add('L', value=L_scaled, min=L_scaled * 1e-1, max=L_scaled * 1e1)

        # If there is no resonance present, we need to determine the value of the inductance via the dataset
        else:
            # Calculate parasitic inductance, scale and write to parameter set
            L = self.calc_L_electrolytic_cap(self.freq, self.z21_data)
            L_scaled = np.mean(L)/config.INDUNIT
            param_set.add('L', value=L_scaled, min=L_scaled*1e-1, max = L_scaled*1e1)

        self.parameters = param_set
        return param_set

    @staticmethod
    def calc_L_electrolytic_cap(freq, data) -> float:
        """
        Function to calculate the parasitic inductance of an electrolytic capacitor

        :param freq: Frequency vector
        :param data: Complex impedance data (needs to be same length as freq)

        :return: The calculated value of the parasitic inductance

        """

        # Calculate omega and the theoretical values for L, assuming that the inductance is present in the entire data
        w = 2*np.pi*freq
        L_test = np.imag(data)/w

        # Reject values smaller than zero (for apparent reasons)
        L_test_valid = L_test[L_test > 0]

        # Now simply take the median, this has shown to perform well on electrolytic caps
        return np.median(L_test_valid)

    def fit_hi_C_model(self, param_set: lmfit.Parameters = None) -> lmfit.Parameters:
        """
        Function to fit the bathtub model

        :param param_set: (optional) A lmfit Parameter set. If not supplied, 'self.parameters' will be used.

        :return: parameters of the fit; also stores the parameters in the instance variable 'self.parameters'

        """

        if param_set is None:
            param_set = self.parameters

        data = self.z21_data
        fit_order = self.order
        mode = constants.fcnmode.FIT_LOG

        modelfreq = self.freq[self.freq<=config.FREQ_UPPER_LIMIT]
        modeldata = data[self.freq<=config.FREQ_UPPER_LIMIT]

        if self.order:
            fit_main_resonance = 0
            out = minimize(self._calculate_Z, param_set,
                           args=(modelfreq, modeldata, fit_order, fit_main_resonance, mode,),
                           method='powell', options={'xtol': 1e-18, 'disp': True})
        else:
            fit_main_resonance = 1
            out = minimize(self._calculate_Z, param_set,
                           args=(modelfreq, modeldata, fit_order, fit_main_resonance, mode,),
                           method='powell', options={'xtol': 1e-18, 'disp': True})

        self.parameters = out.params
        return out.params

    def create_nominal_parameters(self, param_set: lmfit.Parameters = None) -> lmfit.Parameters:
        """
        Function to create parameters for the elements of the main resonance.

        That is L, C, R_Fe and R_s for inductors and C, L, R_iso and R_s for capacitors

        :param param_set: (optional) A lmfit.Parameters() object that the parameters will be written to. If not supplied
            instance variable self.parameters will be used
        :return: A Parameters() object containing the main resonance parameters. Also writes to self.parameters instance
            variable
        """

        if param_set is None:
            param_set = self.parameters

        # Calculate value at the main resonance
        res_value = abs(self.z21_data[self._f0_index])
        w0 = self.f0 * 2 * np.pi

        match self.fit_type:
            case constants.El.INDUCTOR:
                bw_value = res_value / np.sqrt(2)
                f_lower_index = np.flipud(np.argwhere(np.logical_and(self.freq < self.f0, self.data_mag < abs(bw_value))))[0][0]
                f_upper_index = (np.argwhere(np.logical_and(self.freq > self.f0, self.data_mag < abs(bw_value))))[0][0]
                R_Fe = abs(self.z21_data[self._f0_index])
            case constants.El.CAPACITOR:
                R_Iso = constants.R_ISO_VALUE

        match self.fit_type:
            case constants.El.INDUCTOR:
                # Calculate "ideal" capacitance for this resonance
                cap_ideal = 1 / (self.nominal_value * ((self.f0*2*np.pi) ** 2))

                # Scale parameters to units used
                cap_ideal = cap_ideal / config.CAPUNIT
                w0 = w0 / config.FUNIT
                # Create expression string for L
                expression_string_L = '(1/((C*'+ str(config.CAPUNIT)+')*(w0*'+str(config.FUNIT)+')**2))/'+str(config.INDUNIT)

                # Add to parameters
                param_set.add('w0', value=w0, vary=False)
                param_set.add('R_Fe', value=R_Fe, min=constants.MIN_R_FE, max=constants.MAX_R_FE, vary=True)
                param_set.add('R_s', value=self.series_resistance, min=self.series_resistance * 0.1,
                                    max=self.series_resistance * 1.111, vary=True)
                param_set.add('C', value=cap_ideal,
                                    min=cap_ideal * constants.MAIN_RES_PARASITIC_LOWER_BOUND,
                                    max=cap_ideal * constants.MAIN_RES_PARASITIC_UPPER_BOUND, vary=True)
                param_set.add('L', expr = expression_string_L, vary=False)

            case constants.El.CAPACITOR:
                # Calculate ideal inductance for this resonance
                ind_ideal = 1 / (self.nominal_value * ((self.f0 * 2 * np.pi) ** 2))

                # Scale parameters to units used
                ind_ideal = ind_ideal / config.INDUNIT
                w0 = w0 / config.FUNIT
                # Generate expression string for C
                expression_string_C = '(1/((L*' + str(config.INDUNIT) + ')*(w0*' + str(config.FUNIT) + ')**2))/' + str(config.CAPUNIT)

                # Add to parameters
                param_set.add('w0', value=w0, vary=False)
                param_set.add('R_iso', value=R_Iso, min=constants.MIN_R_ISO, max=constants.MAX_R_ISO, vary=True)
                param_set.add('R_s', value=self.series_resistance, min=self.series_resistance * 0.01,
                                    max=self.series_resistance * 1.111, vary=False)
                param_set.add('L', value=ind_ideal,
                                    min=ind_ideal * constants.MAIN_RES_PARASITIC_LOWER_BOUND,
                                    max=ind_ideal * constants.MAIN_RES_PARASITIC_UPPER_BOUND, vary=True)
                param_set.add('C', expr = expression_string_C, vary = False)

        self.parameters = param_set
        return param_set

    def create_higher_order_parameters(self, param_set: lmfit.Parameters = None) -> lmfit.Parameters:
        """
        Method to create the circuit elements for the higher order resonances.
        Can create parameters in two configurations.

        GUI_config 1 constrains one of the energy storing elements by the Q-factor and the other is bound by the resonant
        frequency, leaving only the center frequency and the resistor adjustable.
        GUI_config 2 only constrains one energy storing element to the other by resonant frequency. This is the configuration
        known to perform better if the main resonance fit is fairly accurate.

        :param config_number: selecting the configuration of the resulting parameters. Can be either 1 or 2.
        :param param_set: A Parameters() object that the circuit elements will be written to
        :return: A Parameters() object containing the higher order circuit elements
        """

        ############################## INITIALIZATION ##################################################################

        # Check if custom parameter set has been provided
        if param_set is None:
            param_set = self.parameters

        # Check if the number of detected resonances is higher than the maximum order allowed
        if MAX_ORDER >= len(self.bandwidths):
            self.order = len(self.bandwidths)
        else: # Restrict order to maximum order allowed
            self.order = MAX_ORDER
            self.logger.info("Info: more resonances detected than maximum order permits, "
                             "set order to {value}".format(value=self.order))


        # Initialize array for the modeled bandwidths
        self.modeled_bandwidths = np.zeros([self.order, 3])

        ############################## MAIN LOOP #######################################################################

        # Iterate through all detected resonances; note that we iterate by index since we need a key number
        for key_number in range(1, self.order + 1):

            # Create keys for the parameters to be written
            C_key   = "C%s" % key_number
            L_key   = "L%s" % key_number
            R_key   = "R%s" % key_number
            w_key   = "w%s" % key_number
            BW_key  = "BW%s" % key_number

            # Get upper and lower frequencies of the band
            f_lower = self.bandwidths[key_number - 1][0]
            f_center = self.bandwidths[key_number - 1][1]
            f_upper = self.bandwidths[key_number - 1][2]

            ############################## BANDWIDTH MODEL #############################################################
            # handle bandwidths here -> since the bandwidth detection relies on the 3dB points, which are not always
            # present, we need to model the bandwidth (this is done by brute-force stepping in a separate function)

            # Get indices of the band
            f_center_index = np.where(np.isclose(self.freq, f_center))[0][0]
            f_lower_index = np.where(np.isclose(self.freq, f_lower))[0][0]
            f_upper_index = np.where(np.isclose(self.freq, f_upper))[0][0]

            # Calculate an offset so the bandwidth model receives a bit more datapoints than needed
            # this is done relative to the bandwidth since the datapoints have log-spacing
            n_pts_offset = int(np.floor(((f_upper_index - f_lower_index) / 2) * BW_MODEL_DATA_OFFSET_STRETCH))
            f_lower_index = f_center_index - n_pts_offset
            f_upper_index = f_center_index + n_pts_offset

            # Get data for bandwidth model. Note that we are operating with the smoothed data here
            freq_BW_mdl = self.freq[f_lower_index:f_upper_index]
            data_BW_mdl = self.data_mag[f_lower_index:f_upper_index]*np.exp(1j*np.radians(self.data_ang[f_lower_index:f_upper_index]))

            # invoke bandwidth model
            [f_lower,f_upper,r_value,value_ind,value_cap] = self.model_bandwidth(freq_BW_mdl,data_BW_mdl,f_center)

            #rewrite the obtained bandwidth
            self.modeled_bandwidths[key_number - 1][0] = f_lower
            self.modeled_bandwidths[key_number - 1][1] = f_center
            self.modeled_bandwidths[key_number - 1][2] = f_upper

            # center frequency (omega)
            w_c = f_center * 2 * np.pi
            min_w = w_c * constants.MIN_W_FACTOR
            max_w = w_c * constants.MAX_W_FACTOR

            ############################# BANDWIDTH MODEL PARAMETER CORRECTION #########################################
            # adjust the value of R; since the BW model provides us with an R for the "standalone" circuit, we need
            # to correct it to account for the model data as well, since the model has a non-zero impedance at the
            # point of the newly introduced resonance

            # Calculate impedance data for all resonances we already have, except for the one in question
            curve_data = self._calculate_Z(param_set, self.freq, 2, key_number - 1, 0, constants.fcnmode.OUTPUT)
            data_here = self.z21_data[np.where(np.isclose(self.freq, f_center))[0][0]]
            model_here = curve_data[np.where(np.isclose(self.freq, f_center))[0][0]]
            w_c = f_center * 2 * np.pi
            Q = f_center / (f_upper - f_lower)

            if self.fit_type == constants.El.CAPACITOR:
                R_adjusted = abs(1 / (1 / data_here - 1 / model_here))
                C_adjusted = 1 / (R_adjusted * w_c * Q)

            if self.fit_type == constants.El.INDUCTOR:
                R_adjusted = abs( abs(data_here) -  abs(model_here) )
                C_adjusted = Q / (R_adjusted * w_c)

            r_value = R_adjusted
            value_cap = C_adjusted

            #################### WRITE TO PARAMETER SET ################################################################

            # Rescale parameters to match units
            value_cap = value_cap / config.CAPUNIT
            w_c = w_c / config.FUNIT
            min_w = min_w / config.FUNIT
            max_w = max_w / config.FUNIT

            # Calculate bandwidth value
            #TODO: the bandwidth value should be dependent on the parameters, hence an expression string
            # would be more fitting
            BW_value = (f_upper - f_lower)
            BW_value = BW_value / config.FUNIT

            # Fetch parameter bounds from constants file
            if self.fit_type == El.CAPACITOR:
                # Fetch parameter bounds for capacitors from config
                r_min = r_value * RMINFACTOR_CAP
                r_max = r_value * RMAXFACTOR_CAP
                min_cap = value_cap * CMINFACTOR_CAP
                max_cap = value_cap * CMAXFACTOR_CAP

            elif self.fit_type == El.INDUCTOR:
                # Fetch parameter bounds for inductors from config
                r_min = r_value * RMINFACTOR_COIL
                r_max = r_value * RMAXFACTOR_COIL
                min_cap = value_cap * CMINFACTOR_COIL
                max_cap = value_cap * CMAXFACTOR_COIL

            else:
                #TODO: fit type should be defined at this point in the program; unsure whether or not an exception
                # should be raised here since the edge case should NOT be invoked
                raise AttributeError('Fit type not defined')

            # define expression string for L; L is thereby bound to the corresponding C via w0 = 1/sqrt(L*C)
            # Note that we have to take care of units here; also note that L is dependent on C for all parameters
            expression_string_L = '(1/(('+ C_key +'*'+ str(config.CAPUNIT)+')*('+w_key+'*'+str(config.FUNIT)+')**2))/'+str(config.INDUNIT)

            # finally add parameters with their bounds to param_set
            param_set.add(BW_key, value=BW_value, vary=False)
            param_set.add(C_key, min=min_cap, max=max_cap, value=value_cap, vary=True)
            param_set.add(R_key, value=r_value, min=r_min, max=r_max, vary=True)
            param_set.add(w_key, min=min_w, max=max_w, value=w_c, vary=True)
            param_set.add(L_key, expr=expression_string_L)

        #TODO: write log message for verbose mode
        #TODO: write state of fitter to instance
        self.parameters = param_set
        return param_set

    def pre_fit_bands(self, param_set: lmfit.Parameters = None) -> lmfit.Parameters:
        """
        Method to fit the resonant circuits one-by-one to the impedance data

        The resonant circuits are stepped through. The resonant circuit selected is freed, letting its parameters vary
        while the other parameters are fixed in place. After the fit has completed, the circuits parameters are constrained
        again to not vary and the next circuit will be fit.
        This improves the overall accuracy of the model.

        :param param_set: (optional) A Parameters() object containing the resonant circuits; if not supplied, this
            method takes the instance variable self.parameters
        :return: A Parameters() object containing the fitted bands; also writes the parameters to instance variable
            self.parameters
        """

        # Check if parameters have been supplied, otherwise take instance parameters
        if param_set is None:
            param_set = self.parameters

        # First fix all parameters in place
        self.fix_parameters(param_set)

        for it, band in enumerate(self.modeled_bandwidths):
            # "Cut" the data and frequency vector, so the fitter only looks at the band in question
            fit_freq = self.freq[np.logical_and(self.freq > band[0], self.freq < band[2])]
            fit_data = self.z21_data[np.logical_and(self.freq > band[0], self.freq < band[2])]
            # Generate keys
            key_number = it + 1
            R_key = 'R%s' % key_number
            L_key = 'L%s' % key_number
            C_key = 'C%s' % key_number
            w_key = 'w%s' % key_number

            # Set the parameters in question to vary
            if not param_set[R_key].expr: self.change_parameter(param_set, R_key, vary=True)
            if not param_set[C_key].expr: self.change_parameter(param_set, C_key, vary=True)
            if not param_set[L_key].expr: self.change_parameter(param_set, L_key, vary=True)
            if not param_set[w_key].expr: self.change_parameter(param_set, w_key, vary=True)

            # # Error calculation method for the fit:
            # # choose a log-difference calculation for capacitors and a linear difference calculation for inductors
            # if self.fit_type == El.CAPACITOR:
            #     fit_weighing = fcnmode.FIT_LOG
            # else:
            #     fit_weighing = fcnmode.FIT


            # Do the fit
            out = minimize(self._calculate_Z, param_set,
                           args=(fit_freq, fit_data, self.order, 0, config.FIT_BY,),
                           method='powell', options={'xtol': 1e-18, 'disp': True})

            # Write fit results to parameters and set their 'vary' to False
            R = out.params[R_key].value
            L = out.params[L_key].value
            C = out.params[C_key].value
            w = out.params[w_key].value

            # Write back to parameter set and change the boundaries to the newly obtained estimates
            self.change_parameter(param_set, R_key, value=R, min=R * 0.9, max=R * 1.1, vary=False)
            self.change_parameter(param_set, L_key, value=L, min=L * 0.9, max=L * 1.1, vary=False)
            self.change_parameter(param_set, C_key, value=C, min=C * 0.9, max=C * 1.1, vary=False)
            self.change_parameter(param_set, w_key, value=w, min=w * 0.9, max=w * 1.1, vary=False)

        # Free parameters for further fitting
        self.free_parameters_higher_order(param_set)
        self.parameters = param_set
        return param_set

    def correct_parameters(self, param_set: lmfit.Parameters = None, change_main = False, num_it = 2) -> lmfit.Parameters:
        """
        Method to correct the parameters of the set. Corrects higher order resonances, but can also correct the main
        resonance parameters if needed.

        When the circuit elements are created, their resistance value is roughly the magnitude at resonance minus the
        magnitude of the main resonance at that point. However as new elements are introduced into the model, the peak
        values of the resonances slightly change, since every circuit influences all other circuits. This method aims to
        correct the resistance values and (inherently to keep the same Q-factor for the circuit) the capacitor.

        Note: does not work for the acoustic resonance of MLCCs

        :param param_set: A Parameters() object containing the higher order circuits and main resonance
        :param change_main: Boolean. Selects whether to correct the main resonance parameters
        :param num_it: Number of iterations of the correction (is 2 by default)
        :return: A Parameters() object containing the corrected parameters

        """

        if param_set is None:
            param_set = self.parameters

        params = param_set

        #TODO: method does not really go by DRY paradigm, maybe overhaul

        for it in range(num_it):
            #at the start of each iteration correct main resonance
            curve_data = self._calculate_Z(params, self.freq, 2, self.order, 0, constants.fcnmode.OUTPUT)
            if self.fit_type == constants.El.INDUCTOR and change_main:
                # Get Q value
                #determine upper and lower 3dB-points
                b_l = np.flipud(np.argwhere(np.logical_and(self.freq < self.f0, abs(self.z21_data) < abs(self.z21_data[self.freq == self.f0][0])/(np.sqrt(2)))))[0][0]
                b_u = np.argwhere(np.logical_and(self.freq > self.f0, abs(self.z21_data) < abs(self.z21_data[self.freq == self.f0][0])/(np.sqrt(2)) ))[0][0]
                w0 = (self.f0*2*np.pi)
                Q_main = self.f0 / (self.freq[b_u] - self.freq[b_l])

                #adjust R_Fe by difference from the data
                R_diff = abs(self.z21_data[self.freq <= self.f0][0]) - params['R_Fe'].value
                R_new = params['R_Fe'].value + R_diff
                #adjust C depending on Q and new R_Fe
                C_new = Q_main / (w0*R_new)
                #C_new needs to be rescaled to match units
                C_new = C_new / config.CAPUNIT
                self.change_parameter(params, 'R_Fe', min = R_new *0.8, max = R_new*1.2, value = R_new, vary = False)
                self.change_parameter(params, 'C', min = C_new*0.8, max = C_new*1.2, value= C_new, vary = False)

            elif self.fit_type == constants.El.CAPACITOR and change_main:
                # get Q value
                b_l = np.flipud(np.argwhere(np.logical_and(self.freq < self.f0, abs(self.z21_data) > abs(self.z21_data[self.freq == self.f0][0]) * (np.sqrt(2)))))[0][0]
                b_u = np.argwhere(np.logical_and(self.freq > self.f0, abs(self.z21_data) > abs(self.z21_data[self.freq == self.f0][0]) * (np.sqrt(2))))[0][0]
                w0 = (self.f0 * 2 * np.pi)
                Q_main = self.f0 / (self.freq[b_u] - self.freq[b_l])

                # adjust R_s by difference from the data
                R_diff = abs(self.z21_data[self.freq <= self.f0][0]) - params['R_s'].value
                R_new = params['R_s'].value + R_diff
                # adjust L depending on Q and new R_Fe
                L_new = (R_new*Q_main)/w0
                #rescale L_new to match units
                L_new = L_new / config.INDUNIT
                self.change_parameter(params, 'R_s', min=R_new * 0.8, max=R_new * 1.2, value=R_new, vary=False)
                self.change_parameter(params, 'L', min=L_new * 0.8, max=L_new * 1.2, value=L_new, vary=False)

            for key_number in range(1, self.order + 1):
                index = key_number - 1
                if self.fit_type == constants.El.INDUCTOR:
                    curve_data = self._calculate_Z(params, self.freq, 2, self.order, 0, constants.fcnmode.OUTPUT)
                    band = self.bandwidths[index]
                    dataindex = np.argwhere(self.freq == band[1])[0][0]
                    #check if parameter needs to be corrected -> 5% relative errror is the metric
                    if abs(abs(self.z21_data[dataindex]) - abs(curve_data[dataindex])) / abs(self.z21_data[dataindex]) >= 0.05:

                        R_key = 'R%s' % key_number
                        C_key = 'C%s' % key_number

                        #to get a better estimate for the resistive value, look at the real part and find the peak
                        try:

                            r_peak = abs(self.z21_data[dataindex])

                            # data_lim = self.z21_data[np.logical_and(self.freq < band[2], self.freq > band[0])]
                            # peak = find_peaks(np.real(data_lim), height = 0)
                            # match self.fit_type:
                            #     case constants.El.INDUCTOR:
                            #         r_peak_index = np.argwhere(peak[1]['peak_heights'] == max(peak[1]['peak_heights']))[0][0]
                            #         r_peak = peak[1]['peak_heights'][r_peak_index]
                            #     case constants.El.CAPACITOR:
                            #         r_peak_index = np.argwhere(peak[1]['peak_heights'] == min(peak[1]['peak_heights']))[0][0]
                            #         r_peak = peak[1]['peak_heights'][r_peak_index]
                        except:
                            r_peak = np.real(self.z21_data[dataindex])

                        R = params[R_key].value
                        R_diff = r_peak - abs(curve_data[dataindex])
                        if (R + R_diff) > 0:
                            R_adjusted = R + R_diff
                            w_c = params['w%s' % key_number].value * config.FUNIT
                            BW = params['BW%s' % key_number].value * config.FUNIT
                            Q = w_c / (BW*2*np.pi)
                            C_adjusted = Q / (R_adjusted * w_c)
                            #rescale parameter to match units
                            C_adjusted = C_adjusted / config.CAPUNIT
                            self.change_parameter(params, R_key, min = R_adjusted*0.9, max = R_adjusted*1.1, value = R_adjusted, vary = True, expr ='')
                            self.change_parameter(params, C_key, min = C_adjusted*0.8, max = C_adjusted *1.2, value = C_adjusted, vary = True)
                        else:
                            #if we can't find a valid correction, leave it be
                            pass


        self.parameters = params
        return params

    def _calculate_Z(self, parameters, frequency_vector, data, fit_order, fit_main_res, modeflag):
        """
        Objective function that is invoked by the optimizer. Calculates the impedance for the model.

        This calculates the impedance for the main resonance and then steps through all the higher order resonances in
        order to calculate the impedance for the whole model.

        Has different output modes:

        - Impedance (OUTPUT)
        - Difference in magnitude for the given data (FIT)
        - Difference in phase for the given data (ANGLE)
        - Difference in real part for the given data (FIT_REAL)
        - Difference in imag part for the given data (FIT_IMAG)

        Note: this **does** account for the acoustic resonance of MLCCs.

        :param parameters: A Parameters() object containing the parameters of the model
        :param frequency_vector: The frequency vector over which the impedance is requested.
        :param data: The data to be used for all kinds of FIT output
        :param fit_order: The order of the model (i.e. the number of resonant circuits)
        :param fit_main_res: Boolean. Decides if only the main resonance shall be fit (TRUE) or if higher order
                resonances shall be calculated too (FALSE)
        :param modeflag: Decides which output the function shall give. Can be FIT, OUTPUT, ANGLE, FIT_REAL or FIT_IMAG
        :return: Either the impedance of the model or a difference, depends on modeflag
        """

        #if we only want to fit the main resonant circuit, set order to zero to avoid "for" loops
        if fit_main_res:
            order = 0
        else:
            order = fit_order

        #create array for frequency
        freq = frequency_vector
        w = freq * 2 * np.pi

        #get parameters for main circuit
        C = parameters['C'].value * config.CAPUNIT
        L = parameters['L'].value * config.INDUNIT
        R_s = parameters['R_s'].value

        match self.fit_type:
            case 1:
                R_Fe = parameters['R_Fe'].value
            case 2:
                R_iso = parameters['R_iso'].value

        #calculate main circuits resistance
        XC = 1 / (1j * w * C)
        XL = 1j * w * L
        Z = 0
        match self.fit_type:
            case El.INDUCTOR: #INDUCTOR
                Z_part1 = 1 / ((1 / R_Fe) + (1 / XL))
                Z_main = 1 / ((1 / (R_s + Z_part1)) + (1 / XC))
            case El.CAPACITOR: #CAPACITOR
                Z_main = (1 / ((1 / R_iso) + (1 / XC))) + XL + R_s

        Z = Z_main

        #if MLCC
        if self.captype == constants.captype.MLCC and not fit_main_res:
            Z_A = parameters['R_A'].value + (1j*w*parameters['L_A'].value*config.INDUNIT) + (1/(1j*w*parameters['C_A'].value*config.CAPUNIT))
            Z = 1/(1/Z_main + 1/Z_A)


        for actual in range(1, order + 1):
            key_number = actual
            C_key = "C%s" % key_number
            L_key = "L%s" % key_number
            R_key = "R%s" % key_number
            C_act = parameters[C_key].value * config.CAPUNIT
            L_act = parameters[L_key].value * config.INDUNIT
            R_act = parameters[R_key].value
            Z_C   = 1 / (1j * w * C_act)
            Z_L   = (1j * w * L_act)
            Z_R   = R_act
            match self.fit_type:
                case constants.El.INDUCTOR:
                    Z    += 1 / ( (1/Z_C) + (1/Z_L) + (1/Z_R) )
                case constants.El.CAPACITOR:
                    Z = 1 / ( 1/Z + 1/(Z_R + Z_L + Z_C))

        match modeflag:
            case fcnmode.FIT:
                diff = abs(data)-abs(Z)
                return (diff)
            case fcnmode.FIT_LOG:
                diff = np.log10(abs(data)) - np.log10(abs(Z))
                return (diff)
            case fcnmode.OUTPUT:
                return Z
            case fcnmode.ANGLE:
                return np.angle(data)-np.angle(Z)
            case fcnmode.FIT_REAL:
                return np.real(data)-np.real(Z)
            case fcnmode.FIT_IMAG:
                return np.imag(data)-np.imag(Z)

    def fit_curve_higher_order(self, param_set: lmfit.Parameters = None) -> lmfit.Parameters:
        """
        Method to fit all higher order resonances.

        :param param_set: (optional) A Parameters() object containing the lumped element representation of the model; if
            not supplied, self.parameters is taken as the parameter set
        :return: A Parameters() object containing the fitted circuits; also writes the parameters to instance variable
            self.parameters
        """
        if param_set is None:
            param_set = self.parameters

        # Frequency limit data for fit; only look at the higher order resonances part, offset from main resonance is
        # defined in the constants, maximum frequency is defined in the config
        fit_data_frq_lim = self.z21_data[np.logical_and(self.freq > self.f0 * constants.MIN_ZONE_OFFSET_FACTOR,
                                                   self.freq < config.FREQ_UPPER_LIMIT)]
        freq_data_frq_lim = self.freq[np.logical_and(self.freq > self.f0 * constants.MIN_ZONE_OFFSET_FACTOR,
                                                self.freq < config.FREQ_UPPER_LIMIT)]

        # Fit the parameter set, given that we have an order, otherwise just pass back the parameter set
        if self.order:
            fit_main_resonance = 0
            out = minimize(self._calculate_Z, param_set,
                           args=(freq_data_frq_lim, fit_data_frq_lim, self.order, fit_main_resonance, config.FIT_BY,),
                           method='powell', options={'xtol': 1e-18, 'disp': True})

            self.parameters = out.params
            return out.params
        else:
            self.parameters = param_set
            return param_set

    def fit_main_res_inductor_file_1(self, param_set: lmfit.Parameters = None) -> lmfit.Parameters:
        """
        Method to fit the main resonance circuit for an inductor for the first file (i.e. the reference file)

        :param param_set: A Parameters() object containing the parameters for the main resonance circuit
        :return: A Parameters() object containing the fitted parameters of the main resonance
        """

        if param_set is None:
            param_set = self.parameters


        if constants.DEBUG_FIT: #debug plot -> main res before fit
            self.plot_curve(param_set, 0, 1)

        # Frequency limit data (upper bound) so there are (ideally) no higher order resonances in the main res fit data
        fit_main_resonance = 1
        freq_for_fit = self.freq[(self.freq < self.f0 * constants.MIN_ZONE_OFFSET_FACTOR)]
        data_for_fit = self.z21_data[(self.freq < self.f0 * constants.MIN_ZONE_OFFSET_FACTOR)]

        # Apply offset
        freq_for_fit = freq_for_fit[self._offset:]
        data_for_fit = data_for_fit[self._offset:]

        #################### Main Resonance ############################################################################

        # Start by fitting the main res with all parameters set to vary
        out = minimize(self._calculate_Z, param_set,
                       args=(freq_for_fit, data_for_fit, self.order, fit_main_resonance, config.FIT_BY,),
                       method='powell', options={'xtol': 1e-18, 'disp': True})

        # Set all parameters to not vary; let only R_s vary
        for pname, par in out.params.items():
            par.vary = False
        out.params['R_s'].vary = True

        # Fit R_s via the phase of the data
        out = minimize(self._calculate_Z, out.params,
                       args=(
                           freq_for_fit, data_for_fit, self.order, fit_main_resonance, constants.fcnmode.ANGLE,),
                       method='powell', options={'xtol': 1e-18, 'disp': True})

        # Fitting R_s again does change the main res fit, so set L an C to vary
        out.params['R_s'].vary = False
        out.params['C'].vary = True
        out.params['L'].vary = True

        # And fit again
        out = minimize(self._calculate_Z, out.params,
                       args=(freq_for_fit, data_for_fit, self.order, fit_main_resonance, config.FIT_BY,),
                       method='powell', options={'xtol': 1e-18, 'disp': True})

        # Write series resistance to class variable (important if other files are fit)
        self.series_resistance = out.params['R_s'].value

        # Fix main resonance parameters in place
        self.fix_main_resonance_parameters(out.params)

        # Write to instance variable and return
        self.parameters = out.params
        return out.params

    def fit_main_res_capacitor_file_1(self, param_set: lmfit.Parameters = None) -> lmfit.Parameters:
        """
        Method to fit the main resonance circuit for a capacitor for the first file (i.e. the reference file)

        :param param_set: A Parameters() object containing the parameters for the main resonance circuit
        :return: A Parameters() object containing the fitted parameters of the main resonance
        """

        if param_set is None:
            param_set = self.parameters

        # TODO: think about log-fiting all capacitors
        mode = config.FIT_BY
        if self.captype == constants.captype.HIGH_C:
            mode = constants.fcnmode.FIT_LOG

        fit_main_resonance = 1

        ###################### Main resonance ##########################################################################

        # Frequency limit data (upper bound) so there are (ideally) no higher order resonances in the main res fit data
        freq_for_fit = self.freq[(self.freq < self.f0 * constants.MIN_ZONE_OFFSET_FACTOR)]
        data_for_fit = self.z21_data[(self.freq < self.f0 * constants.MIN_ZONE_OFFSET_FACTOR)]

        # Apply offset
        freq_for_fit = freq_for_fit[self._offset:]
        data_for_fit = data_for_fit[self._offset:]

        # If we have an electrolytic cap, set vary of series resistance to False
        if captype != constants.captype.HIGH_C:
            param_set['R_s'].vary = False

        out = minimize(self._calculate_Z, param_set,
                       args=(freq_for_fit, data_for_fit, self.order, fit_main_resonance, mode,),
                       method='powell', options={'xtol': 1e-18, 'disp': True})

        # Create datasets for data before/after fit
        old_data = self._calculate_Z(param_set, freq_for_fit, [], 0, fit_main_resonance,
                                     constants.fcnmode.OUTPUT)
        new_data = self._calculate_Z(out.params, freq_for_fit, [], 0, fit_main_resonance,
                                     constants.fcnmode.OUTPUT)

        # Check if the fir gives good results, if the fit results are worse than the initial guess take initial guess
        if abs(sum(np.log10(abs(new_data)) - np.log10(abs(data_for_fit)))) < abs(sum(np.log10(abs(old_data)) - np.log10(abs(data_for_fit)))):
            if constants.DEBUG_MESSAGES:
                self.logger.info("CAPFIT MR: did fit")
            param_set = out.params
        else:
            if constants.DEBUG_MESSAGES:
                self.logger.info("CAPFIT MR: did not fit (used param set as is)")

            # This is kind of a hotfix; if there are no higher order resonances present, out.params is not overwritten
            # so the "worse" param config is taken in the final result
            out.params = param_set

        # Fix main resonance parameters in place
        self.fix_main_resonance_parameters(param_set)

        self.parameters = param_set
        return param_set

    #TODO: those two methods are rather redundant (fit_cap/coil_file_n) to the point where they are essentially the
    # same function
    def fit_main_res_inductor_file_n(self, param_set: lmfit.Parameters = None) -> lmfit.Parameters:
        """
        Method to fit the main resonance circuit for an inductor for a file that is not the reference file

        :param param_set: A Parameters() object containing the parameters for the main resonance circuit
        :return: A Parameters() object containing the fitted parameters of the main resonance
        """
        if param_set is None:
            param_set = self.parameters

        mode = config.FIT_BY
        fit_main_resonance = 1

        # Frequency limit data (upper bound) so there are (ideally) no higher order resonances in the main res fit data
        freq_for_fit = self.freq[(self.freq < self.f0 * constants.MIN_ZONE_OFFSET_FACTOR)]
        data_for_fit = self.z21_data[(self.freq < self.f0 * constants.MIN_ZONE_OFFSET_FACTOR)]

        # Apply offset
        freq_for_fit = freq_for_fit[self._offset:]
        data_for_fit = data_for_fit[self._offset:]

        # Fit main resonance
        out = minimize(self._calculate_Z, param_set,
                       args=(freq_for_fit, data_for_fit, self.order, fit_main_resonance, config.FIT_BY,),
                       method='powell', options={'xtol': 1e-18, 'disp': True})

        # Fix main resonance parameters in place
        self.fix_main_resonance_parameters(out.params)

        # Write to instance variable and return
        self.parameters = out.params
        return out.params

    def fit_main_res_capacitor_file_n(self, param_set: lmfit.Parameters = None) -> lmfit.Parameters:
        """
        Method to fit the main resonance circuit for a capacitor for a file that is not the reference file

        :param param_set: A Parameters() object containing the parameters for the main resonance circuit
        :return: A Parameters() object containing the fitted parameters of the main resonance
        """

        if param_set is None:
            param_set = self.parameters

        fit_main_resonance = 1

        # Frequency limit data (upper bound) so there are (ideally) no higher order resonances in the main res fit data
        freq_for_fit = self.freq[(self.freq < self.f0 * constants.MIN_ZONE_OFFSET_FACTOR)]
        data_for_fit = self.z21_data[(self.freq < self.f0 * constants.MIN_ZONE_OFFSET_FACTOR)]

        # Apply offset
        freq_for_fit = freq_for_fit[self._offset:]
        data_for_fit = data_for_fit[self._offset:]

        # Fit main resonance
        out = minimize(self._calculate_Z, param_set,
                       args=(freq_for_fit, data_for_fit, self.order, fit_main_resonance, config.FIT_BY,),
                       method='powell', options={'xtol': 1e-18, 'disp': True})

        # Fix main resonance parameters in place
        self.fix_main_resonance_parameters(out.params)

        # Write to instance variable and return
        self.parameters = out.params
        return param_set

    def select_param_set(self, params: list):
        """
        A method to select the parameter set that fits "better" if multiple parameter sets have been created (that is the
        case if there are multiple configurations)

        This function calls the calculate_band_norm() function which returns a metric of how well a certain parameter set
        fits the bands. That means the parameter sets are tested on how well they fit the corresponding resonances in the
        data.

        :param params: A list of Parameters() objects for the different models to test
        :param debug: Boolean. If True the function logs which set it took
        :return: The selected Parameters() object
        """

        fit_main_resonance = False


        model_data = []
        norm = []
        for it, param_set in enumerate(params):
            model_data.append(self._calculate_Z(param_set, self.freq, [], self.order, fit_main_resonance, constants.fcnmode.OUTPUT))
            norm.append(self.calculate_band_norm(model_data[it]))

        least_norm_mdl_index = np.argwhere(norm == min(norm))[0][0]

        if constants.DEBUG_MESSAGES:
            self.logger.info(self.name + ": selected parameter set " + str(least_norm_mdl_index + 1))

        return params[least_norm_mdl_index]

    def overwrite_main_res_params_file_n(self, param_set0: lmfit.Parameters, param_set: lmfit.Parameters = None)\
            -> lmfit.Parameters:
        """
        Function to overwrite and recalculate the main resonance parameters for a file that is not the reference file.

        :param param_set: A Parameters() object containing the main resonance parameters for the given file
        :param param_set0: A Parameters() object containing the main resonance parameters for the reference file
        :return: A Parameters() object containing the updated main resonance parameters
        """

        if param_set is None:
            param_set = self.parameters

        match self.fit_type:
            case constants.El.INDUCTOR:
                L_ideal = 1 / ((self.f0 * 2 * np.pi) ** 2 * (param_set0['C'].value * config.CAPUNIT))
                self.nominal_value = L_ideal
                L_ideal = L_ideal / config.INDUNIT
                C_val = param_set0['C'].value
                R_s   = param_set0['R_s'].value
                self.change_parameter(param_set, param_name='C', value=C_val, min=C_val * 0.8, max=C_val * 1.2,
                                      vary=False)
                self.change_parameter(param_set, param_name='R_s', value=R_s, min=R_s * 0.8, max=R_s * 1.2,
                                      vary=False)
                self.change_parameter(param_set, param_name='L', value=L_ideal, vary=True, min=L_ideal * 0.8,
                                      max=L_ideal * 1.25, expr='')

            case constants.El.CAPACITOR:
                C_ideal = 1 / ((self.f0 * 2 * np.pi) ** 2 * (param_set0['L'].value * config.INDUNIT))
                self.nominal_value = C_ideal
                C_ideal = C_ideal/config.CAPUNIT
                L_val = param_set0['L'].value
                R_iso = param_set0['R_iso'].value
                # write back value for L and R_iso
                self.change_parameter(param_set, param_name='L', value=L_val, min=L_val * 0.8, max=L_val * 1.2,
                                      vary=False)
                self.change_parameter(param_set, param_name='R_iso', value=R_iso, vary=False)
                # get new value for R_s
                R_s = abs(self.z21_data[self.freq <= self.f0][0])
                self.change_parameter(param_set, param_name='R_s', value=R_s, vary=False)
                self.change_parameter(param_set, param_name='C', value=C_ideal, vary=True, min=C_ideal * 0.999,
                                      max=C_ideal * 1.001, expr='')

        self.parameters = param_set
        return param_set

    def add_higher_order_resonances_MR_fit(self, order: int, param_set0: lmfit.Parameters,
                                           param_set: lmfit.Parameters = None) -> lmfit.Parameters:
        """
        Method to add the higher order parameters of another parameter set to the fitter. This is needed for the main
            resonance fit, in order to have the parameters present for plot generation. It does not serve another
            purpose to add the parameters to the instance variable
        """
        if param_set is None:
            param_set = self.parameters

        self.order = order

        for key_number in range(1, order + 1):
            C_key = "C%s" % key_number
            L_key = "L%s" % key_number
            R_key = "R%s" % key_number
            w_key = "w%s" % key_number
            BW_key = "BW%s" % key_number

            param_set.add(C_key, value=param_set0[C_key].value)
            param_set.add(L_key, value=param_set0[L_key].value)
            param_set.add(R_key, value=param_set0[R_key].value)
            param_set.add(w_key, value=param_set0[w_key].value)
            param_set.add(BW_key, value=param_set0[BW_key].value)

        # Write to instance parameters and return
        self.parameters = param_set
        return param_set



    ####################################V AUXILLIARY V##################################################################
    def write_model_data(self, param_set, model_order):
        """
        Auxilliary function to calculate the impedance of the model

        :param param_set: A Parameters() object containing the circuits of the model
        :param model_order: The order of the model i.e. the number of resonant circuits
        :return: None (stores the model data in an instance variable)
        """
        freq = self.freq
        order = model_order
        fit_main_resonance = 0
        self.model_data = self._calculate_Z(param_set, freq, [], order, fit_main_resonance, constants.fcnmode.OUTPUT)

    @staticmethod
    def change_parameter(param_set, param_name, min=None, max=None, value=None, vary=None, expr=None):
        """
        Auxilliary function to change a single parameter inside a Parameters() object.

        Note: changing one of the attributes to "None" does not work.

        :param param_set: A Parameters() object containing the parameter that will be changed
        :param param_name: The parameter that will be changed
        :param min:  The new min for the parameter. If not given, this attribute will not be changed.
        :param max: The new max for the parameter. If not given, this attribute will not be changed.
        :param value: The new value for the parameter. If not given, this attribute will not be changed.
        :param vary: The new vary state for the parameter. If not given, this attribute will not be changed.
        :param expr: The new expression string for the parameter. If not given, this attribute will not be changed.
        :return: None
        """

        if not min is None:
            param_set[param_name].min = min
        if not max is None:
            param_set[param_name].max = max
        if not value is None:
            param_set[param_name].value = value
        if not vary is None:
            param_set[param_name].vary = vary
        if not expr is None:
            param_set[param_name].expr = expr

    def calculate_band_norm(self, model):
        """
        Auxillary function used in the select_param_set() method. Calculates the difference between the data and the
        model, but only in the region around the resonances.

        :param model: A vector containing the values of the impedance of the model
        :return: cumulative difference (larger cum.diff. means worse fit)
        """
        # Setup
        freq = self.freq
        cumnorm = 0
        zone_factor = 1.2

        # Check the bandwidth regions and check their least squares diff
        for it, band in enumerate(self.bandwidths):
            bandmask = np.logical_and((freq > band[0]/zone_factor), (freq < band[2]*zone_factor))
            raw_data  = abs(self.z21_data[bandmask])
            mdl1_data = abs(model[bandmask])
            norm1 = np.linalg.norm(raw_data - mdl1_data)
            cumnorm += norm1

        return cumnorm

    def plot_curve(self, param_set = None, order = None, main_res = 0, title = None, angle = False):
        """
        Auxillary function to generate a plot of the given measurement data and the impedance of the model.

        :param param_set: A Parameters() object containing the circuits of the model
        :param order: The order of the model (i.e. the number of resonance circuits)
        :param main_res: Boolean. If True the function will only plot the main resonance
        :param title: Optional. String to use as title for the plot
        :return: None
        """
        if param_set is None:
            param_set = self.parameters

        if order is None:
            order = self.order

        testdata = self._calculate_Z(param_set, self.freq, 2, order, main_res, constants.fcnmode.OUTPUT)
        if angle:
            plt.figure()
            plt.semilogx(self.freq, np.rad2deg(np.angle(self.z21_data)))
            plt.semilogx(self.freq, np.rad2deg(np.angle(testdata)))
            if title is not None:
                plt.title(title)
        else:
            plt.figure()
            plt.loglog(self.freq, abs(self.z21_data))
            plt.loglog(self.freq, abs(testdata))
            if title is not None:
                plt.title(title)

    def calc_Z_simple_RLC(self, parameters, freq, data, ser_par, mode):
        """
        Auxillary function that calculates the impedance of a single RLC circuit.

        This method is invoked by the model_bandwidth() funtion and it is an objective function like calculate_Z().

        :param parameters: A Parameters() object, containing R, L and C
        :param freq: The frequency vector over which the impedance is required
        :param data: The corresponding data; needed if FIT mode is selected.
        :param ser_par: Whether the RLC is a serial or parallel resonant circuit (1=serial;2=parallel)
        :param mode: Can be FIT or OUTPUT, selects what will be returned.
        :return: if mode == FIT, will return the difference in magnitude between the data and the model;
                if mode == OUTPUT, will return the impedance of the resonant circuit
        """
        w = np.pi*2*freq
        Z_R = parameters['R'].value
        Z_L = parameters['L'].value * 1j * w
        Z_C = 1 / (parameters['C'].value * 1j * w)

        match ser_par:
            case 1:#serial
                Z = Z_R + Z_L + Z_C
            case 2:#parallel
                Z = 1/(1/Z_R + 1/Z_C + 1/Z_L)

        match mode:
            case constants.fcnmode.FIT:
                # diff = (np.real(data) - np.real(Z)) + (np.imag(data) - np.imag(Z))
                # diff = np.linalg.norm(data-Z)
                diff = abs(data) - abs(Z)
                if constants.DEBUG_BW_MODEL_VERBOSE:
                    test_data = self.calc_Z_simple_RLC(parameters, freq, [], ser_par, 2)
                    plt.loglog(freq,abs(test_data))
                return (diff)
            case constants.fcnmode.OUTPUT:
                return Z

    def model_bandwidth(self, freqdata, data, peakfreq):
        """
        Auxillary function to brute-force fit a resonance circuit around a peak.

        :param freqdata: The frequency vector for said peak
        :param data: The measurement data around the peak
        :param peakfreq: the resonant frequency of the peak
        :return: Tuple: [b_l, b_u, R, L, C] The lower and upper 3dB points from the model (in Hz) and the R, L and C
            values for a circuit that models the band
        """

        #TODO: IMPORTANT RESCALE PARAMETERS!!!!

        #get the height of the peak and the index(will be used later)
        peakindex = np.argwhere(freqdata >= peakfreq)[0][0]
        peakheight = abs(data[peakindex])
        r_val = abs(peakheight)

        #set the flag for parallel/serial circuit
        match self.fit_type:
            case constants.El.INDUCTOR:
                ser_par_flag = 2
            case constants.El.CAPACITOR:
                ser_par_flag = 1

        #find pits in data
        match self.fit_type:
            case constants.El.INDUCTOR:
                pits = find_peaks(-abs(data),prominence=1e-3)
            case constants.El.CAPACITOR:
                pits = find_peaks(abs(data),prominence=1e-3)

        #get the indices of the pits closest to the peak (if they exist)
        try:
            lower_pit_index = pits[0][np.flipud(np.argwhere(pits[0]<peakindex))[0][0]]
        except:
            lower_pit_index = 0

        try:
            upper_pit_index = pits[0][(np.argwhere(pits[0]>peakindex))[0][0]]
        except:
            upper_pit_index = len(data)-1

        #crop data
        modelfreq = freqdata[lower_pit_index:upper_pit_index]
        modeldata = data[lower_pit_index:upper_pit_index]
        #rewrite the peak index after cropping
        peakindex = np.argwhere(modelfreq >= peakfreq)[0][0]

        try:
            #find the inflection points before and behind the peak in order to crop the peak
            # the derivatives need to be sav_gol filtered, otherwise they are too noisy
            ddt = np.gradient(abs(modeldata))
            ddt = scipy.signal.savgol_filter(ddt, 31,3)
            ddt2 = np.gradient(ddt)
            ddt2 = scipy.signal.savgol_filter(ddt2, 31, 3)
            signchange = np.argwhere(abs(np.diff(np.sign(ddt2))))

            #find left and right index of inflection points and crop the data (also constrain the data between 0 and max)
            left = np.clip((signchange[(signchange < peakindex)][-1] - 1), a_min = 0)
            right = np.clip((signchange[(signchange > peakindex)][0] + 1), a_max = len(modelfreq))
            modelfreq = modelfreq[left:right]
            modeldata = modeldata[left:right]
        except:
            #if the derivative approach does not work, we just do it the way it was before (this has less accuracy overall,
            # but should work for most cases)
            pass

        #space through the usual capacitance values in logarithmic steps
        C_max_exp = 3
        C_min_exp = -15
        numsteps = ((C_max_exp-C_min_exp) + 1) + 25
        C_values = np.flipud(np.logspace(C_min_exp, C_max_exp, num=numsteps))
        #C has to be set to some value, so the expr_string for L works
        C = C_values[-1]

        temp_params = Parameters()

        w_c2 = (peakfreq*2*np.pi)**2

        expr_string_L = '1/(C*'+ str(w_c2)+')'
        temp_params.add('R', value=abs(r_val), min = abs(r_val)*0.8,max=abs(r_val)*1.25,vary=False)
        # temp_params.add('L', value = L,  min = L*1e-3, max = L*1e3)
        temp_params.add('C',value = C)#, min = C*1e-3, max = C*1e6)
        temp_params.add('L',expr=expr_string_L)


        ################################################################################################################
        if constants.DEBUG_BW_MODEL_VERBOSE:
            plt.figure()
            plt.loglog(freqdata,abs(data))
            plt.ylim([min(abs(data))-0.5, max(abs(data))+0.5])
        ################################################################################################################


        # now step through the C values and look at the diff from the objective function in order to obtain a good
        # initial guess for lsq fitting

        diff_array = []
        for it,C_val in enumerate(C_values):
            temp_params.add('C',value = C_val, min = C_val*1e-3, max = C_val*1e6)
            diff_data = self.calc_Z_simple_RLC(temp_params, modelfreq, modeldata, ser_par_flag, 1)
            diff_array.append(sum((diff_data)))

        #check if we have a zero crossing
        if any(np.signbit(diff_array) == True):
            sign_change_index = np.argwhere(np.signbit(diff_array) == True)[0][0]
            #do an interpolation in order to find the cap value at the zero crossing -> most accurate value for C
            x = np.linspace(C_values[sign_change_index - 1], C_values[sign_change_index], 10000)
            y = np.linspace(diff_array[sign_change_index - 1], diff_array[sign_change_index], 10000)
            sign_change_index_interp = np.argwhere(np.signbit(y) == True)[0][0]
            C_val_rough_fit = x[sign_change_index_interp]

        else:
            #TODO: handle this case; should not be invoked, when cap values are stepped through well
            # maybe throw an exception and all that... and think about what BW to take then
            pass

        #TODO: look into what to do if the max and min values are too close to each other
        temp_params.add('C',value=C_val_rough_fit, min=C_val_rough_fit * 0.1, max=C_val_rough_fit * 10)

        #do a fit then after we have the approximate value of the cap
        out = minimize(self.calc_Z_simple_RLC, temp_params, args=(modelfreq,modeldata,ser_par_flag,1),
                            method='powell', options={'xtol': 1e-18, 'disp': True})

        ################################################################################################################
        # #PLOTS ( for when you are in the mood for visual analysis ¯\_(ツ)_/¯ )
        if constants.DEBUG_BW_MODEL:
            test_data_again = self.calc_Z_simple_RLC(out.params,freqdata,[],ser_par_flag,2)
            test_data_again_rough = self.calc_Z_simple_RLC(temp_params,freqdata,[],ser_par_flag,2)
            plt.figure()
            plt.plot(diff_array, marker = "D")
            plt.figure()
            plt.loglog(freqdata,abs(data))
            plt.loglog(freqdata,abs(test_data_again))
            plt.loglog(freqdata,abs(test_data_again_rough))

        ################################################################################################################

        #now get the bandwidth
        freq_interp = np.linspace(min(freqdata)-min(freqdata)*(1/1.5),max(freqdata)+max(freqdata)*1.5,num= 10000)
        match self.fit_type:
            case constants.El.INDUCTOR:
                data_interp = self.calc_Z_simple_RLC(out.params, freq_interp, [], ser_par_flag, constants.fcnmode.OUTPUT)
                BW_3_dB_height = peakheight * (1/np.sqrt(2))

                #get the 3dB-Points of the modeled curve
                b_u = freq_interp[np.argwhere(np.logical_and(abs(data_interp) < BW_3_dB_height, freq_interp > peakfreq))[0][0]]
                b_l = freq_interp[np.argwhere(np.logical_and(abs(data_interp) < BW_3_dB_height, freq_interp < peakfreq))[-1][0]]

                return [b_l, b_u, out.params['R'].value, out.params['L'].value, out.params['C'].value]




            case constants.El.CAPACITOR:
                data_interp = self.calc_Z_simple_RLC(out.params, freq_interp, [], ser_par_flag, constants.fcnmode.OUTPUT)

                BW_3_dB_height = peakheight * np.sqrt(2)

                # get the 3dB-Points of the modeled curve
                b_u = freq_interp[np.argwhere(np.logical_and(abs(data_interp) > BW_3_dB_height, freq_interp > peakfreq))[0][0]]
                b_l = freq_interp[np.argwhere(np.logical_and(abs(data_interp) > BW_3_dB_height, freq_interp < peakfreq))[-1][0]]

                return [b_l, b_u, out.params['R'].value, out.params['L'].value, out.params['C'].value]

    def fix_main_resonance_parameters(self, param_set):
        """
        Auxillary function to lock the main resonance parameters in place

        :param param_set: A Parameters() object containing the main resonance parameters that shell be fixed
        :return: None
        """
        param_set['R_s'].vary = False
        param_set['L'].vary = False
        param_set['C'].vary = False
        match self.fit_type:
            case constants.El.INDUCTOR:
                param_set['R_Fe'].vary = False
            case constants.El.CAPACITOR:
                param_set['R_iso'].vary = False

    def free_parameters_higher_order(self, param_set):
        """
        Auxillary function to free (i.e. set vary to True) the higher order circuit parameters.

        :param param_set: A Parameters object, containing the higher order resonances
        :return: None
        """

        for key_number in range(1, self.order + 1):
            #create keys
            C_key   = "C%s" % key_number
            L_key   = "L%s" % key_number
            R_key   = "R%s" % key_number
            w_key   = "w%s" % key_number
            param_set[C_key].vary = True
            param_set[L_key].vary = True
            param_set[R_key].vary = True
            param_set[w_key].vary = True

    def fix_parameters(self, param_set, R = True, L=True, C=True, w=True):
        #TODO: make this a private function maybe?
        """
        Auxilliary function to lock parameters in a Parameters() object

        :param param_set: A Parameters() object containing the parameters to lock
        :param R: Boolean. Whether to lock these parameters or not. True locks the parameters. Will lock all higher order Rs
        :param L: Boolean. Whether to lock these parameters or not. True locks the parameters. Will lock all higher order Ls
        :param C: Boolean. Whether to lock these parameters or not. True locks the parameters. Will lock all higher order Cs
        :param w: Boolean. Whether to lock these parameters or not. True locks the parameters. Will lock all higher order ws
        :return: None
        """

        param_set['R_s'].vary = False
        param_set['C'].vary = False
        param_set['L'].vary = False
        match self.fit_type:
            case constants.El.CAPACITOR:
                param_set['R_iso'].vary = False
            case constants.El.INDUCTOR:
                param_set['R_Fe'].vary = False
        for key_number in range(1, self.order + 1):
            #create keys
            C_key   = "C%s" % key_number
            L_key   = "L%s" % key_number
            R_key   = "R%s" % key_number
            w_key   = "w%s" % key_number
            param_set[C_key].vary = not C
            param_set[L_key].vary = not L
            param_set[R_key].vary = not R
            param_set[w_key].vary = not w

    ############################## LOGGING AND WARNING #################################################################
    def _debug_info(self, message: str):
        """
        Function to print debug info

        Parameters
        ----------
        message

        Returns
        -------

        """
        pass

    def _warning(self, message:str):
        pass

    def _log_info(self, message: str):
        """
        Function to print a log message

        Parameters
        ----------
        message

        Returns
        -------

        """
        pass

    ################################V EXPERIMENTAL V####################################################################


