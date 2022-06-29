
# The fitter class shall take the task of fitting the data, as well as smoothing it and performing manipulations
# This is likely to become a rather long task, especially for CMCs and this class is therefore likely to be long
# I do not know yet what it will have to contain and how to best handle the data
# Most of this class will be based on Payer's program
#NOTE: THE CLASS IN THE FORM THAT IT IS NOW IS NOT ABLE TO MANAGE MULTIPLE FILES!!!!!
import copy

import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy import signal as sg
from lmfit import minimize, Parameters
from scipy.signal import find_peaks
import decimal
import sys
import pandas as pd


# import constants
import fitterconstants
from fitterconstants import *


class Fitter:

    def __init__(self, logger_instance):
        self.nominal_value = None
        self.parasitive_resistance = None
        self.prominence = None
        self.saturation = None
        self.file = None
        self.z21_data = None
        self.data_mag = None
        self.data_ang = None
        self.model_data = None

        self.fit_type = None
        self.ser_shunt = None

        self.out = None

        self.parameters = Parameters()

        self.logger = logger_instance

        self.frequency_zones = None
        self.bandwidths = None
        self.bad_bandwidth_flag = None
        self.peak_heights = None
        self.frequency_vector = None



        self.f0 = None
        self.max_order = fitterconstants.MAX_ORDER
        self.order = 0



    ####################################################################################################################
    # Parsing Methods
    ####################################################################################################################

    #method to set the entry values of the specification
    def set_specification(self, pass_val, para_r, prom, sat, fit_type):

        self.fit_type = fit_type

        if pass_val is None:
            #if we do not have the nominal value try to calculate it
            try:
                self.calculate_nominal_value()
            #if we can't calculate it, pass the exception back to the calling function
            except Exception as e:
                raise e
        else:
            self.nominal_value = pass_val

        if para_r is None:
            self.calculate_nominal_Rs()
        else:
            self.parasitive_resistance = para_r

        if prom is None:
            self.prominence = fitterconstants.PROMINENCE_DEFAULT
        else:
            self.prominence = prom

        if sat is None:
            self.saturation = None
        else:
            self.saturation = sat



    #method to parse the files from the iohandler
    def set_file(self, file):
        self.file = file
        try:
            self.frequency_vector = self.file.data.f
            self.logger.info("File: " + self.file.name)
        except Exception:
            raise Exception("No Files were provided, please select a file!")




    ####################################################################################################################
    # Pre-Processing Methods
    ####################################################################################################################

    def calc_series_thru(self, Z0):
        self.z21_data = 2 * Z0 * ((1 - self.file.data.s[:, 1, 0]) / self.file.data.s[:, 1, 0])
        self.ser_shunt = fitterconstants.calc_method.SERIES

    def calc_shunt_thru(self, Z0):
        self.z21_data = (Z0 * self.file.data.s[:, 1, 0]) / (2 * (1 - self.file.data.s[:, 1, 0]))
        self.ser_shunt = fitterconstants.calc_method.SHUNT

    def crop_data(self,crop):
        self.z21_data = self.z21_data[crop:]
        self.frequency_vector = self.frequency_vector[crop:]

    def smooth_data(self):
        # Use Savitzky-Golay filter for smoothing the input data, because in the region of the global minimum there is
        # oscillation. After filtering a global minimum can be found easier.
        sav_gol_mode = 'interp'
        self.data_mag = scipy.signal.savgol_filter(abs(self.z21_data), fitterconstants.SAVGOL_WIN_LENGTH,
                                                   fitterconstants.SAVGOL_POL_ORDER, mode=sav_gol_mode)
        self.data_ang = scipy.signal.savgol_filter(np.angle(self.z21_data, deg=True), fitterconstants.SAVGOL_WIN_LENGTH,
                                                   fitterconstants.SAVGOL_POL_ORDER, mode=sav_gol_mode)
        #limit the data to +/- 90°
        self.data_ang = np.clip(self.data_ang, -90, 90)

        return 0

    def calculate_nominal_value(self):
        offset = fitterconstants.NOMINAL_VALUE_CALC_OFFSET
        nominal_value = 0
        freq = self.frequency_vector

        match self.fit_type:
            case El.INDUCTOR:
                # find first point where the phase crosses 0 using numpy.argwhere --> f0
                index_angle_smaller_zero = np.argwhere(self.data_ang < 0)
                index_ang_zero_crossing = index_angle_smaller_zero[0][0]
                f0 = freq[index_ang_zero_crossing]

                # if the first zero crossing is smaller that the offset i.e. the first zero crossing is at the start of
                # the data, raise an exception
                if index_ang_zero_crossing <= offset:
                    raise Exception("Error: Could not calculate nominal value;"
                                    "the Phase of the dataset seems to be bad, consider cropping the data")

                if max(self.data_ang[offset:index_ang_zero_crossing]) < 85:
                    #if we can't detect the nominal value raise exception
                    raise Exception("Error: Inductive range not detected (max phase = {value}°).\n"
                                    "Please specify nominal inductance.".format(value=np.round(max(self.data_ang), 1)))

                #crop data to [offset:f0] in order to find the linear range for the calculation of nominal value
                ang_test_data = self.data_ang[freq < f0][offset:]
                curve_data = self.z21_data[freq < f0][offset:]
                w_data = (freq[freq < f0][offset:])*2*np.pi

                #then ceil the data, and select all values that have max value (we can assume linear range here, ideally
                # phase should be 90° here)
                # ang_test_ceil = np.ceil(ang_test_data)
                # bool_select = 1#ang_test_ceil == max(ang_test_ceil)

                #create an array filled with possible values for L; calculation is L = imag(Z)/w
                L_vals = []
                for it, curve_sample in enumerate(zip(curve_data, w_data)):
                    #if bool_select[it]:
                    L_vals.append(np.imag(curve_sample[0])/curve_sample[1])

                #find the 50% quantile of the slope data and define the max slope allowed
                quantile_50 = np.quantile(np.gradient(self.data_mag)[freq<f0],0.5)
                max_slope = quantile_50 * fitterconstants.QUANTILE_MULTIPLICATION_FACTOR
                #boolean index the data that has lower than max slope and calculate the mean
                L_vals_eff = np.array(L_vals)[np.gradient(self.data_mag)[freq<f0][offset:] < max_slope]
                self.nominal_value = np.mean(L_vals_eff)

                output_dec = decimal.Decimal("{value:.3E}".format(value=self.nominal_value)) #TODO: this has to be normalized output to 1e-3/-6/-9 etc
                self.logger.info("Nominal Inductance not provided, calculated: " + output_dec.to_eng_string())


            case El.CAPACITOR:
                # find first point where the phase crosses 0
                index_angle_larger_zero = np.argwhere(self.data_ang > 0)
                index_ang_zero_crossing = index_angle_larger_zero[0][0]
                f0 = freq[index_ang_zero_crossing]

                # if the first zero crossing is smaller that the offset i.e. the first zero crossing is at the start of
                # the data, raise an exception
                if index_ang_zero_crossing <= offset:
                    raise Exception("Error: Could not calculate nominal value;"
                                    "the Phase of the dataset seems to be bad, consider cropping the data")

                if min(self.data_ang[offset:index_ang_zero_crossing]) > -85:
                    raise Exception("Error: Capacitive range not detected (min phase = {value}°).\n"
                                    "Please specify nominal capacitance.".format(value=np.round(min(self.data_ang), 1)))

                #TODO: the calculation for nominal parameters for capacitors seems to be off --> check on that again

                #crop data to [offset:f0] in order to find the linear range for the calculation of nominal value
                ang_test_data = self.data_ang[freq < f0][offset:]
                curve_data = self.z21_data[freq < f0][offset:]
                w_data = (freq[freq < f0][offset:])*2*np.pi

                #then floor the data, and select all values that have max value (we can assume linear range here, ideally
                # phase should be -90° here)
                # ang_test_ceil = np.floor(ang_test_data)
                # bool_select = ang_test_ceil == min(ang_test_ceil)

                #create an array filled with possible values for L; calculation is L = imag(Z)/w
                C_vals = []
                for it, curve_sample in enumerate(zip(curve_data, w_data)):
                    # if bool_select[it]:
                    C_vals.append(-1/(np.imag(curve_sample[0])*curve_sample[1]))


                # find the 50% quantile of the slope data and define the max slope allowed
                quantile_50 = np.quantile(np.gradient(self.data_mag)[freq < f0], 0.5)
                max_slope = quantile_50 * fitterconstants.QUANTILE_MULTIPLICATION_FACTOR
                # boolean index the data that has lower than max slope and calculate the mean
                C_vals_eff = np.array(C_vals)[np.gradient(self.data_mag)[freq < f0][offset:] < max_slope]
                self.nominal_value = np.mean(C_vals_eff)



                # test_values = []
                # for sample in range(offset, len(freq)):
                #     if self.data_ang[sample] == min(self.data_ang[offset:index_ang_zero_crossing]):
                #         nominal_value = 1 / (2 * np.pi * freq[sample] * self.data_mag[sample])
                #         self.nominal_value = nominal_value
                #         test_values.append(self.nominal_value)
                #         # break
                # try:
                #     test_values_gradient = abs(np.gradient(test_values, 2))
                #     # it takes the first values instead of the "linear" range. need to fix this. possibly by taking the longest min gradient
                #     #TODO: i dont have a single clue how this works for capacitors EDIT: if we can't calculate gradient, just don't do it i guess
                #     self.nominal_value = test_values[np.argmin(np.amin(test_values_gradient))]
                # except Exception:
                #     pass
                output_dec = decimal.Decimal("{value:.3E}".format(value=self.nominal_value))
                self.logger.info("Nominal Capacitance not provided, calculated: " + output_dec.to_eng_string())

            case 3:
                self.nominal_value = 0



        return self.nominal_value

    def calculate_nominal_Rs(self):
        R_s_input = min(abs(self.z21_data))
        self.parasitive_resistance = R_s_input
        #log
        output_dec = decimal.Decimal("{value:.3E}".format(value=R_s_input))
        self.logger.info("Nominal Resistance not provided, calculated: " + output_dec.to_eng_string())

    def get_main_resonance(self):
        #TODO: this method goes by the phase, it could use some more 'robustness'

        freq = self.frequency_vector

        #set w0 to 0 in order to have feedback, if the method didn't work
        w0 = 0

        match self.fit_type:

            case 1: #INDUCTOR
                index_angle_smaller_zero = np.argwhere(self.data_ang < 0)
                index_ang_zero_crossing = index_angle_smaller_zero[0][0]
                continuity_check = index_angle_smaller_zero[10][0]

            case 2: #CAPACITOR
                index_angle_larger_zero = np.argwhere(self.data_ang > 0)
                index_ang_zero_crossing = index_angle_larger_zero[0][0]
                continuity_check = index_angle_larger_zero[10][0]

            case 3: #CMC
                sign = 1

        # TODO: there could be some problems here: a) the resonant frequency could be at the start of the data and
        #   b) the resonant frequency could be at the end of data... those are cases in which the phase data is faulty.
        #   an exception should be raised already in this case, but only if the calculate nominal value method was run

        if continuity_check:
            f0 = freq[index_ang_zero_crossing]
            w0 = f0 * 2 * np.pi
            self.f0 = f0
            #log and print
            output_dec = decimal.Decimal("{value:.3E}".format(value=f0))
            self.logger.info("Detected f0: "+ output_dec.to_eng_string())
            print("Detected f0: "+output_dec.to_eng_string())

        if w0 == 0:
            raise Exception('ERROR: Main resonant frequency could not be determined.')

    def get_resonances(self): #TODO: tidy up this whole method :/

        R_s = self.parasitive_resistance
        freq = self.frequency_vector
        #create one figure for the resonance plots
        if fitterconstants.DEBUG_BW_DETECTION:
            plt.figure()

        # in order to use the same methods for capacitors as we do for inductors, we simply flip the dataset, so we still
        # detect "peaks" although there are pits
        match self.fit_type:
            case fitterconstants.El.INDUCTOR:
                magnitude_data = self.data_mag
                phase_data = self.data_ang
                peak_min_height = np.log10(R_s)
            case fitterconstants.El.CAPACITOR:
                magnitude_data = self.data_mag
                phase_data = self.data_ang
                #TODO: maybe look into that; there might be a better option for the peak height
                peak_min_height = min(magnitude_data)

        #frequency limit the data
        magnitude_data = magnitude_data[freq < fitterconstants.FREQ_UPPER_LIMIT]
        phase_data = phase_data[freq < fitterconstants.FREQ_UPPER_LIMIT]
        freq = freq[freq < fitterconstants.FREQ_UPPER_LIMIT]


        prominence_mag = self.prominence
        prominence_phase = self.prominence

        #find peaks of Magnitude Impedance curve (using scipy.signal.find_peaks)
        match self.ser_shunt:
            case fitterconstants.calc_method.SERIES:
                mag_maxima = find_peaks(20*np.log10(magnitude_data), height=peak_min_height, prominence=prominence_mag)
            case fitterconstants.calc_method.SHUNT:
                mag_maxima = find_peaks(-20 * np.log10(magnitude_data), height=peak_min_height, prominence=prominence_mag)


        mag_minima = find_peaks(magnitude_data * -1, prominence=prominence_mag)
        #find peaks of Phase curve
        phase_maxima = find_peaks(phase_data, prominence=prominence_phase)
        phase_minima = find_peaks(phase_data * -1, prominence=prominence_phase)

        ######FOR MANUAL TESTING!!!VVV
        # test_prom = find_peaks(magnitude_data, height=peak_min_height, prominence=prominence_mag)
        # plt.loglog(magnitude_data)
        # plt.plot(test_prom[0], test_prom[1]['peak_heights'], marker='D', linestyle='')

        #map to frequency; TODO: we are using the file here, so if there are multiple files, need to change this
        #TODO: why are we even calculating the magnitude maxima if they are never used???
        f_mag_maxima = freq[mag_maxima[0]]
        f_mag_minima = freq[mag_minima[0]]

        f_phase_maxima = freq[phase_maxima[0]]
        f_phase_minima = freq[phase_minima[0]]

        #ignore all peaks that lie "before" the main resonance and that are to close to the main resonance
        min_zone_start = self.f0 * fitterconstants.MIN_ZONE_OFFSET_FACTOR

        # TODO: delete unnecessary variables here
        ang_minima_pos = f_phase_minima[f_phase_minima > min_zone_start]
        ang_maxima_pos = f_phase_maxima[f_phase_maxima > min_zone_start]

        ang_minima_pos = f_mag_minima[f_mag_minima > min_zone_start]
        ang_maxima_pos = f_mag_maxima[f_mag_maxima > min_zone_start]

        mag_minima_pos = f_mag_minima[f_mag_minima > min_zone_start]
        mag_maxima_pos = f_mag_maxima[f_mag_maxima > min_zone_start]

        mag_minima_index = mag_minima[0][f_mag_minima > min_zone_start]
        mag_maxima_index = mag_maxima[0][f_mag_maxima > min_zone_start]

        mag_maxima_value = mag_maxima[1]['peak_heights'][f_mag_maxima > min_zone_start]


        # plot commands to check peak values TODO: this is for testing
        # markerson = mag_maxima[0]
        # plt.loglog(self.data_mag,'-bD', markevery=markerson)
        # plt.show()
        # plt.figure()

        number_zones = len(mag_maxima_pos)
        bandwidth_list = []
        peak_heights = []
        bad_BW_flag = np.zeros((number_zones,2))
        for num_maximum in range(0, number_zones):
            #resonance frequency, corresponding height and index
            res_fq = ang_maxima_pos[num_maximum]
            res_index = mag_maxima_index[num_maximum]
            res_value = magnitude_data[res_index]

            #get 3dB value
            # if we work with an inductor the curve is mirrored, so we have to account for that
            match self.fit_type:
                case fitterconstants.El.INDUCTOR:
                    bw_value = res_value / np.sqrt(2)
                case fitterconstants.El.CAPACITOR:
                    bw_value = res_value*np.sqrt(2)

            try:
                #find the index where the 3db value is reached; also check if the frequency is lower than the resonance,
                #but higher than the min zone; if that does not work use the default offset
                #NOTE: since we need the first value in front of the resonance we have to flipud the array
                f_lower_index = np.flipud(np.argwhere(np.logical_and(freq > min_zone_start, np.logical_and(freq < res_fq, (magnitude_data) < (bw_value)))))[0][0]
            except IndexError:
                f_lower_index = res_index - fitterconstants.DEFAULT_OFFSET_PEAK
                bad_BW_flag[num_maximum][0] = 1

            try:
                f_upper_index = np.argwhere(np.logical_and(freq > res_fq, (magnitude_data) < (bw_value)))[0][0]
            except IndexError:
                #here we need to account for the fact that we could overshoot the max index
                if res_index + fitterconstants.DEFAULT_OFFSET_PEAK < len(freq):
                    f_upper_index = res_index + fitterconstants.DEFAULT_OFFSET_PEAK
                    bad_BW_flag[num_maximum][1] = 1
                else:
                    f_upper_index = len(freq) - 1

            # check if the found 3dB points are in an acceptable range i.e. not "behind" the next peak or "in front of"
            # the previous peak. If that is the case we set the index to a default offset to get a "bandwidth"
            if num_maximum != 0:
                if f_lower_index < mag_maxima_index[num_maximum - 1]:
                    f_lower_index = res_index - fitterconstants.DEFAULT_OFFSET_PEAK
                    bad_BW_flag[num_maximum] = 1
            if num_maximum < number_zones-1:
                if f_upper_index > mag_maxima_index[num_maximum + 1]:
                    #again we could overshoot the max index here
                    if res_index + fitterconstants.DEFAULT_OFFSET_PEAK < len(freq):
                        f_upper_index = res_index + fitterconstants.DEFAULT_OFFSET_PEAK
                        bad_BW_flag[num_maximum][1] = 1
                    else:
                        f_upper_index = len(freq) - 1
                        bad_BW_flag[num_maximum][1] = 1

            # this checks if the value of the upper/lower bound is greater than the value of the resonance peak
            # that is the case if we chose the default offset #TODO: look into how to handle this case

            # if ((magnitude_data[res_index]) < (magnitude_data[f_upper_index])) or ((magnitude_data[res_index]) < (magnitude_data[f_lower_index])):
            #     # at the moment we are just skipping the peak in that case
            #     pass
            # else:
            f_tuple = [freq[f_lower_index], res_fq, freq[f_upper_index]]
            bandwidth_list.append(f_tuple)
            peak_heights.append(abs(res_value))
            #THIS IS FOR TESTING
            if fitterconstants.DEBUG_BW_DETECTION:
                markerson = [f_lower_index,res_index,f_upper_index]
                plt.loglog(self.data_mag, '-bD', markevery=markerson)

        try:
            #spread BW of last circuit; TODO: maybe center the band?
            stretch_factor = fitterconstants.BANDWIDTH_STRETCH_LAST_ZONE
            # bandwidth_list[-1][0] = max(freq) * (1/strech_factor)

            # bandwidth_list[-1][2]= max(freq)*5

            #
            bandwidth_list[-1][2] = bandwidth_list[-1][2] * stretch_factor
            bandwidth_list[-1][0] = bandwidth_list[-1][0] * stretch_factor

            #peak_heights[-1] = abs(max(self.z21_data)) * 2
        except IndexError:
            self.logger.info("INFO: No resonances found except the main resonance, consider a lower value for the prominence")


        #TODO: this block is for testing the bandwidth model

        # # testing the modeled bandwidth here
        # mdl_offset = 0
        # for it in range(0, number_zones):
        #     if bad_BW_flag[it]:
        #         freq_data = freq[np.logical_and(freq > bandwidth_list[it][0],freq < bandwidth_list[it][2])]
        #         mdl_mag_data = magnitude_data[np.logical_and(freq > bandwidth_list[it][0],freq < bandwidth_list[it][2])]
        #         mdl_phase_data = phase_data[np.logical_and(freq > bandwidth_list[it][0], freq < bandwidth_list[it][2])]
        #         mdl_data = mdl_mag_data * np.exp(1j * np.radians(mdl_phase_data))
        #         self.model_bandwidth(freq_data,mdl_data)

        self.peak_heights = peak_heights
        self.bandwidths = bandwidth_list
        self.bad_bandwidth_flag = bad_BW_flag

    def create_nominal_parameters(self):

        self.parameters.add('R_s', value=self.parasitive_resistance, min=self.parasitive_resistance*0.9,max=self.parasitive_resistance*1.111, vary=False)




        #get bandwidth

        freq = self.frequency_vector
        res_value = self.z21_data[freq == self.f0]
        w0 = self.f0 * 2 * np.pi
        self.parameters.add('w0',value=w0,vary=False)

        match self.fit_type:
            case fitterconstants.El.INDUCTOR:
                bw_value = res_value / np.sqrt(2)
                f_lower_index = np.flipud(np.argwhere(np.logical_and(freq < self.f0, self.data_mag < bw_value)))[0][0]
                f_upper_index = (np.argwhere(np.logical_and(freq > self.f0, self.data_mag < bw_value)))[0][0]
                BW = freq[f_upper_index] - freq[f_lower_index]
                R_Fe = (self.f0 * (self.f0 * 2 * np.pi) * self.nominal_value) / BW
            case fitterconstants.El.CAPACITOR:
                bw_value = res_value * np.sqrt(2)
                f_lower_index = np.flipud(np.argwhere(np.logical_and(freq < self.f0, self.data_mag > bw_value)))[0][0]
                f_upper_index = (np.argwhere(np.logical_and(freq > self.f0, self.data_mag > bw_value)))[0][0]
                BW = freq[f_upper_index] - freq[f_lower_index]
                R_Iso = fitterconstants.R_ISO_VALUE#BW/(self.f0 * (self.f0 * 2 * np.pi)*self.nominal_value)





        match self.fit_type:
            case fitterconstants.El.INDUCTOR:
                #calculate "perfect" capacitor for this resonance
                cap_ideal = 1 / (self.nominal_value * ((self.f0*2*np.pi) ** 2))
                #add to parameters

                expression_string_L = '1/(' + str(self.f0*2*np.pi) + '**2*' + 'C)'

                self.parameters.add('R_Fe', value=R_Fe, min=fitterconstants.MIN_R_FE, max=fitterconstants.MAX_R_FE, vary=True)

                # #config A
                # #main element
                # self.parameters.add('L', value=self.nominal_value, min=self.nominal_value * 0.9, max=self.nominal_value * 1.1, vary=False)
                # # expression_string_C = '1/( w0 **2* L)'
                # self.parameters.add('C', value=cap_ideal,
                #                     min=cap_ideal * fitterconstants.MAIN_RES_PARASITIC_LOWER_BOUND,
                #                     max=cap_ideal * fitterconstants.MAIN_RES_PARASITIC_UPPER_BOUND, vary=True)

                #Config B
                self.parameters.add('C', value=cap_ideal,
                                    min=cap_ideal * fitterconstants.MAIN_RES_PARASITIC_LOWER_BOUND,
                                    max=cap_ideal * fitterconstants.MAIN_RES_PARASITIC_UPPER_BOUND, vary=True)
                self.parameters.add('L', expr = expression_string_L, vary=False)


                #alternative -> varies the main element and keeps the parasitic element constrained via expression
                # self.parameters.add('L', value=self.nominal_value, min=self.nominal_value * 0.9,max=self.nominal_value * 1.1, vary=True)
                # self.parameters.add('C',expr=expression_string_C, vary=False)


            case fitterconstants.El.CAPACITOR:
                # calculate "perfect" inductor for this resonance
                ind_ideal = 1 / (self.nominal_value * ((self.f0 * 2 * np.pi) ** 2))

                self.parameters.add('R_iso', value=R_Iso, min=fitterconstants.MIN_R_ISO, max=fitterconstants.MAX_R_ISO, vary=True)
                #main element


                expression_string_C = '1/(' + str(self.f0*2*np.pi) + '**2*' + 'L)'
                self.parameters.add('L', value=ind_ideal,
                                    min=ind_ideal * fitterconstants.MAIN_RES_PARASITIC_LOWER_BOUND,
                                    max=ind_ideal * fitterconstants.MAIN_RES_PARASITIC_UPPER_BOUND, vary=True)

                self.parameters.add('C', expr = expression_string_C, vary = False)

                # self.parameters.add('L',expr=expression_string_L,vary=False)


            case 3:
                #TODO: CMCs -> eh scho wissen
                dummy = 0

        return 0

    def create_elements(self, config_number):

        #if we got too many frequency zones -> restrict fit to max order
        #else get order from frequency zones and write found order to class
        if self.max_order >= len(self.bandwidths):
            order = len(self.bandwidths)
            self.order = len(self.bandwidths)
        else:
            order = self.max_order
            self.order = order
            self.logger.info("Info: more resonances detected than maximum order permits, set order to {value}".format(value=order))
            #TODO: and also throw and except please
            #TODO: some methods are not robust enough for this fit maybe?


        #model the main res, so we have good estimates for the impedance curve
        #NOTE: THIS DOES NOT WORK THAT EASY VVVVVVVVVVVVV

        # f_c_index = np.argwhere(self.frequency_vector == self.f0)[0][0]
        # f_l_index = 0
        # f_u_index = np.argwhere(self.frequency_vector == self.f0)[0][0]+1000
        # # get data for bandwidth model
        # freq_BW_mdl = self.frequency_vector[f_l_index:f_u_index]
        # data_BW_mdl = self.data_mag[f_l_index:f_u_index] * np.exp(1j * np.radians(self.data_ang[f_l_index:f_u_index]))
        # # now model the BW
        # [b_l, b_u, r_value, value_ind, value_cap] = self.model_bandwidth(freq_BW_mdl, data_BW_mdl, self.f0)

        main_res_data = self.calculate_Z(self.parameters, self.frequency_vector,2,0,1,fitterconstants.fcnmode.OUTPUT)


        for key_number in range(1, order + 1):

            #create keys
            C_key   = "C%s" % key_number
            L_key   = "L%s" % key_number
            R_key   = "R%s" % key_number
            w_key   = "w%s" % key_number
            BW_key  = "BW%s" % key_number




            #get upper and lower frequencies
            b_l = self.bandwidths[key_number - 1][0]
            b_c = self.bandwidths[key_number - 1][1]
            b_u = self.bandwidths[key_number - 1][2]

            # handle bandwidths here -> since the bandwidth detection relies on the 3dB points, which are not always
            # present, we may need to "model" the BW. If we have one of the two 3dB points though, we can assume symmetric
            # bandwidth EDIT: bandwidth model has been applied to all peaks, since it gives good estimates for the
            # parameter values

            stretch_factor = 1.5
            #get indices of the band
            f_c_index = np.argwhere(self.frequency_vector == b_c)[0][0]
            f_l_index = np.argwhere(self.frequency_vector == b_l)[0][0]
            f_u_index = np.argwhere(self.frequency_vector == b_u)[0][0]
            #calculate diffference between upper and lower, so the number of points is relative to where we are in
            #the data, since the measurement points are not equally spaced
            n_pts_offset = ((f_u_index - f_l_index) / 2) * stretch_factor
            #recalc lower and upper bound
            f_l_index = f_c_index - int(np.floor(n_pts_offset))
            f_u_index = f_c_index + int(np.floor(n_pts_offset))
            #get data for bandwidth model
            freq_BW_mdl = self.frequency_vector[f_l_index:f_u_index]
            data_BW_mdl = self.data_mag[f_l_index:f_u_index]*np.exp(1j*np.radians(self.data_ang[f_l_index:f_u_index]))

            #upper and lower 3dB point faulty
            if self.bad_bandwidth_flag[key_number-1].all:
                #now model the BW
                [b_l,b_u,r_value,value_ind,value_cap] = self.model_bandwidth(freq_BW_mdl,data_BW_mdl,b_c)
            #only lower 3dB point faulty
            elif self.bad_bandwidth_flag[key_number-1][0]:
                [_,_, r_value, value_ind, value_cap] = self.model_bandwidth(freq_BW_mdl, data_BW_mdl, b_c)
                b_l = b_c - (b_u - b_c)
            # only upper 3dB point faulty
            elif self.bad_bandwidth_flag[key_number - 1][1]:
                [_,_, r_value, value_ind, value_cap] = self.model_bandwidth(freq_BW_mdl, data_BW_mdl, b_c)
                b_u = b_c + (b_c - b_l)
            #both points present
            else:
                [_,_, r_value, value_ind, value_cap] = self.model_bandwidth(freq_BW_mdl, data_BW_mdl, b_c)


            # bandwidth
            BW_min = (b_u - b_l) * fitterconstants.BW_MIN_FACTOR
            BW_max = (b_u - b_l) * fitterconstants.BW_MAX_FACTOR
            BW_value = (b_u - b_l)  # BW_max / 8

            # calculate Q-factor
            q = b_c / BW_value

            # center frequency (omega)
            w_c = b_c * 2 * np.pi
            min_w = w_c * fitterconstants.MIN_W_FACTOR
            max_w = w_c * fitterconstants.MAX_W_FACTOR

            #TODO: look into how to handle the min==max errors

            # #we get an error of min==max for the parameters, if one parameter is too small-> just subtract something i guess
            # while value_cap < fitterconstants.MINIMUM_PRECISION:
            #     value_cap = value_cap * 2
            #



            r_max = r_value * 1.25
            r_min = r_value * 0.5


            min_ind = value_ind * 2
            max_ind = value_ind *0.5

            #add parameters

            #this is the default configuration, i.e. the config how tristan had it
            # self.parameters.add(BW_key, min=BW_min,     max=BW_max,     value=BW_value              , vary=True)
            # self.parameters.add(w_key,  min=min_w,      max=max_w,      value=w_c                   , vary=True)
            # self.parameters.add(C_key,  min=min_cap,    max=max_cap,    value=value_cap             , vary=True)
            # self.parameters.add(L_key,  min=1e-20,      max=L,          expr=expression_string_L    , vary=False)
            # self.parameters.add(R_key,  min=1e-3,       max=1e4,        expr=expression_string_R    , vary=False)




            if self.fit_type == fitterconstants.El.CAPACITOR:
                # good values for capacitor fitting
                max_cap = value_cap * 1e2
                min_cap = value_cap * 1e-2

                r_max = r_value * 1.01
                r_min = r_value * 0.990

                expression_string_L = '1/(' + w_key + '**2*' + C_key + ')'
                expression_string_R = '(1/(' + w_key + '/(' + BW_key + '*' + str(2 * np.pi) + ')))*sqrt(' + L_key + '/' + C_key + ')'

                self.parameters.add(w_key, min=min_w, max=max_w, value=w_c, vary=False)
                self.parameters.add(BW_key, min=BW_min, max=BW_max, value=BW_value, vary=False)
                self.parameters.add(C_key, min=min_cap, max=max_cap, value=value_cap, vary=True)

                self.parameters.add(L_key, expr=expression_string_L, vary=False)
                self.parameters.add(R_key, expr=expression_string_R, vary=False)

                match config_number:
                    case 1:
                        # config B -> default config; this goes via the Q factor
                        expression_string_L = '1/(' + w_key + '**2*' + C_key + ')'
                        expression_string_R = '((' + BW_key + '*' + str(2 * np.pi) + ')/(' + w_key  + '))*sqrt(' + L_key + '/' + C_key + ')'
                        self.parameters.add(w_key, min=min_w, max=max_w, value=w_c, vary=True)
                        self.parameters.add(L_key, expr=expression_string_L, vary=False)
                        self.parameters.add(R_key, expr=expression_string_R, vary=False)
                    case 2:
                        # config D (assuming perfectly fitted main resonance)
                        if (r_value - abs(main_res_data[f_c_index])) > 0:
                            r_value = r_value - abs(main_res_data[f_c_index]) * (
                                        abs(main_res_data[f_c_index]) / r_value)
                            value_cap = q / (w_c * r_value)
                            max_cap = value_cap * 1e1  # 2
                            min_cap = value_cap * 1e-1  # 500e-3

                        expression_string_L = '1/(' + w_key + '**2*' + C_key + ')'
                        expression_string_C = L_key + '*(' + '(' + w_key + '/(' + BW_key + '*' + str(
                            2 * np.pi) + '))' + '/' + R_key + ')**2'
                        expression_string_C = '((' + w_key + '/(' + BW_key + '*' + str(
                            2 * np.pi) + ')))/(' + R_key + '*' + w_key + ')'
                        self.parameters.add(C_key, min=min_cap, max=max_cap, value=value_cap, vary=True)
                        self.parameters.add(R_key, value=r_value, min=r_value * 0.2, max=r_value * 5)
                        self.parameters.add(w_key, min=min_w, max=max_w, value=w_c, vary=True)
                        self.parameters.add(L_key, expr=expression_string_L)
                        # self.parameters.add(C_key, expr=expression_string_C)






            else:

                max_cap = value_cap * 1e2#2
                min_cap = value_cap * 1e-2#500e-3


                self.parameters.add(w_key, min=min_w, max=max_w, value=w_c, vary=False)

                self.parameters.add(BW_key, min=BW_min, max=BW_max, value=BW_value, vary=False)

                self.parameters.add(C_key, min=min_cap, max=max_cap, value=value_cap, vary=True)

                # #config A -> does not perform too well
                # self.parameters.add(R_key, value=r_value, min=r_min, max=r_max, vary=True)
                # self.parameters.add(L_key, expr=expression_string_L, vary=False)

                match config_number:
                    case 1:
                        #config B -> default config; this goes via the Q factor
                        expression_string_L = '1/(' + w_key + '**2*' + C_key + ')'
                        expression_string_R = '(' + w_key + '/(' + BW_key + '*' + str(2*np.pi) + '))*sqrt(' +  L_key + '/' + C_key + ')'
                        self.parameters.add(w_key, min=min_w, max=max_w, value=w_c, vary=True)
                        self.parameters.add(L_key, expr=expression_string_L, vary=False)
                        self.parameters.add(R_key, expr=expression_string_R, vary=False)
                    case 2:
                        # config D (assuming perfectly fitted main resonance)
                        if (r_value - abs(main_res_data[f_c_index])) > 0:
                            r_value = r_value - abs(main_res_data[f_c_index])*(abs(main_res_data[f_c_index])/r_value)
                            value_cap = q/(w_c*r_value)
                            max_cap = value_cap * 1e1  # 2
                            min_cap = value_cap * 1e-1  # 500e-3

                        expression_string_L = '1/(' + w_key + '**2*' + C_key + ')'
                        expression_string_C = L_key + '*(' + '(' + w_key + '/(' + BW_key + '*' + str(2*np.pi) + '))' + '/' + R_key + ')**2'
                        expression_string_C = '((' + w_key + '/(' + BW_key + '*' + str(2*np.pi) + ')))/('+R_key+'*'+w_key+')'
                        self.parameters.add(C_key, min=min_cap, max=max_cap, value=value_cap, vary=True)
                        self.parameters.add(R_key, value = r_value, min = r_value * 0.2, max = r_value * 5)
                        self.parameters.add(w_key, min=min_w, max=max_w, value=w_c, vary=True)
                        self.parameters.add(L_key, expr=expression_string_L)
                        # self.parameters.add(C_key, expr=expression_string_C)


                # # config C
                # expression_string_L = '((' + R_key + '**2)*' + C_key + ')/(' + str(q ** 2) + ')'
                # self.parameters.add(R_key, value=r_value, min=r_min, max=r_max, vary=True)
                # self.parameters.add(L_key, expr=expression_string_L, vary=False)


        return 0

    def calculate_Z(self, parameters, frequency_vector, data, fit_order, fit_main_res, modeflag):
        #method to calculate the impedance curve from chained parallel resonance circuits
        #this method is needed for the fitter

        #if we only want to fit the main resonant circuit, set order to zero to avoid "for" loops
        if fit_main_res:
            order = 0
        else:
            order = fit_order

        #create array for frequency
        freq = frequency_vector
        w = freq * 2 * np.pi

        #get parameters for main circuit
        C = parameters['C'].value
        L = parameters['L'].value
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

                #trying a different model here
                # Z_main = 1 / ( (1 / R_iso) + (1 / (XC + R_s + XL)) )

        Z = Z_main

        for actual in range(1, order + 1):
            key_number = actual
            C_key = "C%s" % key_number
            L_key = "L%s" % key_number
            R_key = "R%s" % key_number
            C_act = parameters[C_key].value
            L_act = parameters[L_key].value
            R_act = parameters[R_key].value
            Z_C   = 1 / (1j * w * C_act)
            Z_L   = (1j * w * L_act)
            Z_R   = R_act
            match self.fit_type:
                case fitterconstants.El.INDUCTOR:
                    Z    += 1 / ( (1/Z_C) + (1/Z_L) + (1/Z_R) )
                case fitterconstants.El.CAPACITOR:
                    Z = 1 / ( 1/Z + 1/(Z_R + Z_L + Z_C))

        # diff = (np.real(data) - np.real(Z)) + 1j * (np.imag(data) - np.imag(Z))
        # return abs(diff)

        match modeflag:
            case fcnmode.FIT:
                diff = abs(data)-abs(Z)#(np.real(data) - np.real(Z)) + (np.imag(data) - np.imag(Z))
                # diff = (np.real(data) - np.real(Z))**2 + (np.imag(data) - np.imag(Z))**2
                # diff = np.linalg.norm(data-Z)
                return (diff)
            case fcnmode.OUTPUT:
                return Z

    def start_fit_file_1(self):

        freq = self.frequency_vector
        fit_data = self.z21_data
        fit_order = self.order
        mode = fitterconstants.fcnmode.FIT

        if fitterconstants.DEBUG_FIT: #debug plot -> main res before fit
            self.plot_curve_before_fit()
        
        #frequency limit data (upper bound) so there are (ideally) no higher order resonances in the main res fit data
        fit_main_resonance = 1
        freq_for_fit = freq[(freq < self.f0 * fitterconstants.MIN_ZONE_OFFSET_FACTOR)]
        data_for_fit = fit_data[(freq < self.f0 * fitterconstants.MIN_ZONE_OFFSET_FACTOR)]

        #crop some samples of the start of data (~100) because the slope at the start of the dataset might be off
        freq_for_fit = freq_for_fit[fitterconstants.MAIN_RES_FIT_OFFSET_SAMPLES:]
        data_for_fit = data_for_fit[fitterconstants.MAIN_RES_FIT_OFFSET_SAMPLES:]

        #now do the fit
        self.out = minimize(self.calculate_Z, self.parameters,
                            args=(freq_for_fit, data_for_fit, fit_order, fit_main_resonance, mode,),
                            method='powell', options={'xtol': 1e-18, 'disp': True})

        #create datasets for data before/after fit
        old_data = self.calculate_Z(self.parameters, freq, [], 0, fit_main_resonance, fitterconstants.fcnmode.OUTPUT)
        new_data = self.calculate_Z(self.out.params, freq, [], 0, fit_main_resonance, fitterconstants.fcnmode.OUTPUT)

        #TODO: fix fit behaviour--> the worse (before fit) parameter set is taken for the 0v_cap file at the moment
        #check if the main resonance fit yields good results -> else: go with initial guess
        if np.linalg.norm(new_data - self.z21_data) < np.linalg.norm(old_data - self.z21_data):
            self.parameters = self.out.params
        else:
            #redundant, but for readability
            self.parameters = self.parameters

        #fix main resonance parameters in place
        self.fix_main_resonance_parameters()

        if fitterconstants.DEBUG_FIT:#debug plot -> fitted main resonance
            self.plot_curve_before_fit()

        ###################### Higher order resonances #################################################################

        # Frequency limit for fit data
        fit_data_frq_lim = fit_data[np.logical_and(freq > self.f0 * fitterconstants.MIN_ZONE_OFFSET_FACTOR ,freq < fitterconstants.FREQ_UPPER_LIMIT)]
        freq_data_frq_lim = freq[np.logical_and(freq > self.f0 * fitterconstants.MIN_ZONE_OFFSET_FACTOR ,freq < fitterconstants.FREQ_UPPER_LIMIT)]

        ##################### config 1 fit

        #now create the elements for higher order circuits, fit with config 1
        self.create_elements(1)
        fit_order = self.order

        if fitterconstants.DEBUG_FIT:#debug plot -> curve with initial guesses for the elements
            self.plot_curve_before_fit()


        #check if we have higher order resonances, otherwise the fitter can't do anything here (yields error: too many fcn calls)
        if fit_order:
            fit_main_resonance = 0
            self.out = minimize(self.calculate_Z, self.parameters,
                                args=(freq_data_frq_lim, fit_data_frq_lim, fit_order, fit_main_resonance, mode,),
                                method='powell', options={'xtol': 1e-18, 'disp': True})

        param_set_1 = self.out.params

        ##################### config 2 fit

        # now create the elements for higher order circuits, fit with config 1
        self.create_elements(2)
        fit_order = self.order

        if fitterconstants.DEBUG_FIT:  # debug plot -> curve with initial guesses for the elements
            self.plot_curve_before_fit()

        # check if we have higher order resonances, otherwise the fitter can't do anything here (yields error: too many fcn calls)
        if fit_order:
            fit_main_resonance = 0
            self.out = minimize(self.calculate_Z, self.parameters,
                                args=(freq_data_frq_lim, fit_data_frq_lim, fit_order, fit_main_resonance, mode,),
                                method='powell', options={'xtol': 1e-18, 'disp': True})

        param_set_2 = self.out.params

        #calculate model data for both parameter sets
        mode = fitterconstants.fcnmode.OUTPUT
        model_data_1 = self.calculate_Z(param_set_1, freq, [], self.order, fit_main_resonance, mode)
        model_data_2 = self.calculate_Z(param_set_2, freq, [], self.order, fit_main_resonance, mode)

        if fitterconstants.DEBUG_FIT:  # debug plot -> both parameter sets after fit
            plt.figure()
            plt.loglog(self.frequency_vector, abs(fit_data))
            plt.loglog(freq, abs(model_data_1))
            plt.loglog(freq, abs(model_data_2))

        norm_mdl1 = self.calculate_band_norm(model_data_1)
        norm_mdl2 = self.calculate_band_norm(model_data_2)

        #decide which param set to take
        if norm_mdl1 < norm_mdl2:
        # if abs(sum(abs(fit_data_frq_lim)-abs(norm_1_mdl))) < abs(sum(abs(fit_data_frq_lim)-abs(norm_2_mdl))):
            self.parameters = param_set_1
            self.model_data = model_data_1
            if fitterconstants.DEBUG_FIT:
                self.logger.info("Debug: Took Parameter config 1 (Q constrained)")
        else:
            self.parameters = param_set_2
            self.model_data = model_data_2
            if fitterconstants.DEBUG_FIT:
                self.logger.info("Debug: Took Parameter config 2 (R free/C free)")

        ################################################################################################################

        self.parameters.pretty_print()




        #for testin purposes
        if fitterconstants.DEBUG_FIT:
            plt.figure()
            plt.loglog(self.frequency_vector, abs(fit_data))
            plt.loglog(freq, abs(self.model_data))


        return 0

    def start_fit_file_n(self, fitting_mode):
        #fix parameters in place, so the high order resonances are not affected by the fitting process of the current/
        # voltage dependent main element

        freq = self.frequency_vector
        fit_data = self.z21_data
        res_value = self.z21_data[freq == self.f0]

        #determine wether to perform a full fit (i.e. fit all parameters) or if only the main resonance should be fit
        match fitting_mode:
            case fitterconstants.multiple_fit.MAIN_RES_FIT:
                self.fix_parameters()
                fit_main_resonance = 1
            case fitterconstants.multiple_fit.FULL_FIT:
                #if we want to perform a full fit, we unfortunately have to remake the parameters, without loosing the
                #value for C though (R_s is stored as instance variable anyways)
                match self.fit_type:
                    case fitterconstants.El.INDUCTOR:
                        C_val = self.parameters['C'].value
                        R_Fe_val = self.parameters['R_Fe'].value
                        #clear and re-initiate parameters
                        self.parameters = Parameters()
                        self.create_nominal_parameters()
                        self.get_resonances()
                        self.create_elements()
                        #write back value for C and keep it in place
                        self.parameters['C'].value = C_val
                        self.parameters['C'].vary = False
                    case fitterconstants.El.CAPACITOR:
                        L_val = self.parameters['L'].value
                        # clear and re-initiate parameters
                        self.parameters = Parameters()
                        self.create_nominal_parameters()
                        self.get_resonances()
                        self.create_elements()
                        # write back value for C and keep it in place
                        self.parameters['L'].value = L_val
                        self.parameters['L'].vary = False
                fit_main_resonance = 0



        #calculate ideal value for the dependent element, so we are as close as possible to the detected resonance
        match self.fit_type:
            case fitterconstants.El.INDUCTOR:
                L_ideal = 1 / ((self.f0 * 2 * np.pi)**2 * self.parameters['C'].value)
                self.parameters['L'].value = L_ideal
                self.parameters['L'].vary = True
                self.parameters['L'].min = L_ideal*0.8
                self.parameters['L'].max = L_ideal*1.25
                self.parameters['L'].expr = ''
                #
                bw_value = res_value / np.sqrt(2)
                f_lower_index = np.flipud(np.argwhere(np.logical_and(freq < self.f0, self.data_mag < bw_value)))[0][0]
                f_upper_index = (np.argwhere(np.logical_and(freq > self.f0, self.data_mag < bw_value)))[0][0]
                BW = freq[f_upper_index] - freq[f_lower_index]
                R_Fe = (self.f0 * (self.f0 * 2 * np.pi) * self.nominal_value) / BW

                match fitting_mode:
                    case fitterconstants.multiple_fit.FULL_FIT:
                        self.parameters['R_Fe'].vary = True
                        self.parameters['R_Fe'].value = R_Fe
                        self.parameters['R_Fe'].min = R_Fe * 0.8
                        self.parameters['R_Fe'].max = R_Fe * 1.25
                        freq_for_fit = freq
                        data_for_fit = fit_data
                    case fitterconstants.multiple_fit.MAIN_RES_FIT:
                        self.parameters['R_Fe'].vary = False
                        freq_for_fit = freq[(freq < self.f0 * fitterconstants.MIN_ZONE_OFFSET_FACTOR)]
                        data_for_fit = fit_data[(freq < self.f0 * fitterconstants.MIN_ZONE_OFFSET_FACTOR)]
                        freq_for_fit = freq_for_fit[fitterconstants.MAIN_RES_FIT_OFFSET_SAMPLES:]
                        data_for_fit = data_for_fit[fitterconstants.MAIN_RES_FIT_OFFSET_SAMPLES:]


                #
                # self.parameters['R_s'].vary = True

            case fitterconstants.El.CAPACITOR:
                C_ideal = 1 / ((self.f0 * 2 * np.pi) * self.parameters['L'].value)
                self.parameters['C'].value = C_ideal
                self.parameters['C'].vary = True
                self.parameters['C'].min = C_ideal * 0.8
                self.parameters['C'].max = C_ideal * 1.25



        mode = fitterconstants.fcnmode.FIT
        # if only main res fit -> order = 0; fit_main_res = 1
        fit_order = self.order
        # call the minimizer and pass the arguments
        self.out = minimize(self.calculate_Z, self.parameters,
                            args=(freq_for_fit, data_for_fit, fit_order, fit_main_resonance, mode,),
                            method='powell', options={'xtol': 1e-18, 'disp': True})

        self.logger.info("debug: calculated main element{value:.3E}".format(value = self.parameters['L'].value))

        model_data = self.calculate_Z(self.out.params, freq,2,self.order,0,fitterconstants.fcnmode.OUTPUT)
        if fitterconstants.DEBUG_FIT:
            plt.figure()
            plt.loglog(self.frequency_vector, abs(fit_data))
            plt.loglog(freq, abs(model_data))
            # manager = plt.get_current_fig_manager()
            # manager.full_screen_toggle()
        self.out.params.pretty_print()
        # plt.show()


    ####################################V AUXILLIARY V##################################################################


    def calculate_band_norm(self, model):
        #function to compare the two fit results
        freq = self.frequency_vector
        cumnorm = 0
        zone_factor = 1.2

        #check the bandwidth regions and check their least squares diff
        for it, band in enumerate(self.bandwidths):
            bandmask = np.logical_and((freq > band[0]/zone_factor), (freq < band[2]*zone_factor))
            raw_data  = abs(self.z21_data[bandmask])
            mdl1_data = abs(model[bandmask])
            norm1 = np.linalg.norm(raw_data - mdl1_data)
            cumnorm += norm1

        return cumnorm


    def plot_curve_before_fit(self):

        testdata = self.calculate_Z(self.parameters, self.frequency_vector, 2, self.order, 0, 2)
        plt.figure()
        plt.loglog(self.frequency_vector, abs(self.z21_data))
        plt.loglog(self.frequency_vector, abs(testdata))


    def calc_Z_simple_RLC(self,parameters,freq,data,ser_par,mode):
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
            case fitterconstants.fcnmode.FIT:
                # diff = (np.real(data) - np.real(Z)) + (np.imag(data) - np.imag(Z))
                # diff = np.linalg.norm(data-Z)
                diff = abs(data) - abs(Z)
                if fitterconstants.DEBUG_BW_MODEL_VERBOSE:
                    test_data = self.calc_Z_simple_RLC(parameters, freq, [], ser_par, 2)
                    plt.loglog(freq,abs(test_data))
                return (diff)
            case fitterconstants.fcnmode.OUTPUT:
                return Z


    def model_bandwidth(self, freqdata, data, peakfreq):

        #TODO: adapt for capacitors!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        #get the height of the peak and the index(will be used later)
        peakindex = np.argwhere(freqdata == peakfreq)[0][0]
        peakheight = abs(data[peakindex])
        r_val = abs(peakheight)

        #set the flag for parallel/serial circuit
        match self.fit_type:
            case fitterconstants.El.INDUCTOR:
                ser_par_flag = 2
            case fitterconstants.El.CAPACITOR:
                ser_par_flag = 1

        #find pits in data
        match self.fit_type:
            case fitterconstants.El.INDUCTOR:
                pits = find_peaks(-abs(data))
            case fitterconstants.El.CAPACITOR:
                pits = find_peaks(abs(data))

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
        temp_params.add('C',value = C, min = C*1e-3, max = C*1e6)
        temp_params.add('L',expr=expr_string_L)


        ################################################################################################################
        if fitterconstants.DEBUG_BW_MODEL_VERBOSE:
            plt.figure()
            plt.loglog(freqdata,abs(data))
            plt.ylim([min(abs(data))-0.5, max(abs(data))+0.5])
        ################################################################################################################


        # now step through the C values and look at the diff from the objective function in order to obtain a good
        # initial guess for lsq fitting

        diff_array = []
        for C_val in C_values:
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
        if fitterconstants.DEBUG_BW_MODEL:
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
            case fitterconstants.El.INDUCTOR:
                data_interp = self.calc_Z_simple_RLC(out.params, freq_interp, [], ser_par_flag,fitterconstants.fcnmode.OUTPUT)
                BW_3_dB_height = peakheight * (1/np.sqrt(2))

                #get the 3dB-Points of the modeled curve
                b_u = freq_interp[np.argwhere(np.logical_and(abs(data_interp) < BW_3_dB_height, freq_interp > peakfreq))[0][0]]
                b_l = freq_interp[np.argwhere(np.logical_and(abs(data_interp) < BW_3_dB_height, freq_interp < peakfreq))[-1][0]]

                return [b_l, b_u, out.params['R'].value, out.params['L'].value, out.params['C'].value]




            case fitterconstants.El.CAPACITOR:
                data_interp = self.calc_Z_simple_RLC(out.params, freq_interp, [], ser_par_flag, fitterconstants.fcnmode.OUTPUT)

                BW_3_dB_height = peakheight * np.sqrt(2)

                # get the 3dB-Points of the modeled curve
                b_u = freq_interp[np.argwhere(np.logical_and(abs(data_interp) > BW_3_dB_height, freq_interp > peakfreq))[0][0]]
                b_l = freq_interp[np.argwhere(np.logical_and(abs(data_interp) > BW_3_dB_height, freq_interp < peakfreq))[-1][0]]

                return [b_l, b_u, out.params['R'].value, out.params['L'].value, out.params['C'].value]


    def fix_main_resonance_parameters(self):
        self.parameters['R_s'].vary = False
        self.parameters['L'].vary = False
        self.parameters['C'].vary = False
        match self.fit_type:
            case fitterconstants.El.INDUCTOR:
                self.parameters['R_Fe'].vary = False
            case fitterconstants.El.CAPACITOR:
                self.parameters['R_iso'].vary = False


    def fix_parameters(self):
        #method to fix the parameters in place, except the nominal value, which varies with different current/voltage
        match self.fit_type:
            case fitterconstants.El.INDUCTOR:
                self.parameters['C'].vary = False
                self.parameters['R_Fe'].vary = False
                self.parameters['R_s'].vary = False
            case fitterconstants.El.CAPACITOR:
                self.parameters['L'].vary = False
                self.parameters['R_iso'].vary = False
                self.parameters['R_s'].vary = False

        for key_number in range(1, self.order + 1):
            #create keys
            C_key   = "C%s" % key_number
            L_key   = "L%s" % key_number
            R_key   = "R%s" % key_number
            w_key   = "w%s" % key_number
            BW_key  = "BW%s" % key_number
            self.parameters[C_key].vary = False
            self.parameters[L_key].vary = False
            self.parameters[R_key].vary = False
            self.parameters[w_key].vary = False
            self.parameters[BW_key].vary = False


    def overwrite_main_resonance_parameters(self):
        # method to overwrite the nominal parameters with the parameters obtained by modeling the main resonance circuit

        self.parameters['R_s'].value = self.out.params['R_s'].value
        self.parameters['R_s'].vary = False

        match self.fit_type:
            case 1:
                self.parameters['L'].value = self.out.params['L'].value
                self.parameters['R_Fe'].value = self.out.params['R_Fe'].value
                self.parameters['C'].value = self.out.params['C'].value
                self.parameters['L'].vary = False
                self.parameters['R_Fe'].vary = False
                self.parameters['C'].vary = False
            case 2:
                self.parameters['L'].value = self.out.params['L'].value
                self.parameters['R_iso'].value = self.out.params['R_iso'].value
                self.parameters['C'].value = self.out.params['C'].value
                self.parameters['L'].vary = False
                self.parameters['R_iso'].vary = False
                self.parameters['C'].vary = False
            case 3:
                # TODO: CMCs -> eh scho wissen
                dummy = 0

        return 0

    #################################### V OBSOLETE V###################################################################


    def fit_end_zone(self):
        #method to fit the slope at the end of data and fix it... this is just a test, so it might not work as intended
        temp_params = Parameters()
        key_number = self.order
        C_key = "C%s" % key_number
        L_key = "L%s" % key_number
        R_key = "R%s" % key_number
        temp_params.add("C")
        temp_params.add("L")
        temp_params.add("R_Fe")
        temp_params.add("R_s", value = 0, vary = False)

        offset = 2700

        temp_params["C"] = self.parameters[C_key]
        temp_params["L"] = self.parameters[L_key]
        temp_params["R_Fe"] = self.parameters[R_key]

        freq = self.frequency_vector[offset:-1]
        fit_data = self.z21_data[offset:-1]
        fit_order = 0
        fit_main_resonance = 1
        mode = fitterconstants.fcnmode.FIT

        end_zone_model = minimize(self.calculate_Z, temp_params,
                       args=(freq, fit_data, fit_order, fit_main_resonance, mode,),
                       method='leastsq')

        mode = fitterconstants.fcnmode.OUTPUT
        freq = np.linspace(min(self.frequency_vector),max(self.frequency_vector)+1e12,120000)
        model_data = self.calculate_Z(end_zone_model.params,freq, [2], self.order, fit_main_resonance, mode)

        plt.figure()
        plt.loglog(freq, abs(model_data))
        plt.loglog(self.frequency_vector, abs(self.z21_data))
        plt.show()




    def test_fit_main_res(self):
        #function for testing purposes NOT IN USE

        C = self.out.params['C'].value
        L = self.out.params['L'].value
        R_Fe = self.out.params['R_Fe'].value
        R_s = self.out.params['R_s'].value

        w = self.frequency_vector * 2 * np.pi
        XC = 1 / (1j * w * C)
        XL = 1j * w * L
        Z = 0
        Z_part1 = 1 / ((1 / R_Fe) + (1 / XL))
        Z_main = 1 / ((1 / (R_s + Z_part1)) + (1 / XC))
        plt.loglog(Z_main)
        plt.loglog(self.data_mag)


    def fit_iteration_callback(self, out_model):
        #method to set the parameters of the output to the parameters used for the fit TODO:rewrite this comment
        #TODO: delete -> UNUSED METHOD
        for key in self.parameters.keys():
            self.parameters[key].value = out_model.params[key].value
            self.parameters[key].min = self.parameters[key] * 0.8
            self.parameters[key].max = self.parameters[key] * 1.2


    def calculate_transfer_function_data_diff(self,params,freq,data):

        #initiate arrays for numerator and denominator
        a_list = []
        b_list = []
        #calculate omega
        w = freq *2*np.pi
        #get the values from the parameters
        b_list.append(params['b0'].value)
        for it in range(1,self.order):
            a_list.append(params['a%s' % it].value)
            b_list.append(params['b%s' % it].value)
        #create a transfer function and get magnitude and phase
        TF = sg.TransferFunction(a_list,b_list)
        w, mag, ph = TF.bode(w, n=len(self.z21_data))
        #calculate (r+ji)
        TF_data = mag * np.exp(np.radians(ph))

        diff = (np.real(data) - np.real(TF_data)) + (np.imag(data) - np.imag(TF_data))
        return diff


    def fit_transfer_function(self):

        #TODO: this needs to be invoked otherwise the self.order is not set

        self.create_elements()

        params = Parameters()

        a_list = []
        b_list = []

        params.add('b0', value=1, min=1e-6, max=1e6)
        for it in range(1,self.order):
            params.add('a%s' % it, value = 1, min = 1e-6, max = 1e6)
            params.add('b%s' % it, value = 1, min = 1e-6, max = 1e6)
            a_list.append(params['a%s' % it].value)
            b_list.append(params['b%s' % it].value)

        out = minimize(self.calculate_transfer_function_data_diff, params,
                            args=(self.frequency_vector, self.z21_data),
                            method='powell', options={'xtol': 1e-18, 'disp': True})

        #calculate transfer function data
        a_list=[]
        b_list=[]
        b_list.append(out.params['b0'].value)
        for it in range(1,self.order):
            a_list.append(out.params['a%s' % it].value)
            b_list.append(out.params['b%s' % it].value)

        TF = sg.TransferFunction(a_list,b_list)

        w, mag, ph = TF.bode(2*np.pi*self.frequency_vector, n = len(self.z21_data))

        plt.figure()
        plt.loglog(w,self.z21_data)
        plt.loglog(w,mag)



        pass





