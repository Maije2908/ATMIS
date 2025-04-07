import config
import skrf as rf
import pandas as pd
import os
from fitter import *
import constants
import matplotlib
from matplotlib import pyplot as plt

class IOhandler:
    """
    The IOHandler class takes care of the filehandling.
    It loads files and stores a reference to it.
    It also takes care of the output, generating the LTSpice Netlists as well as generating the Bode-plots for output
    """

    def __init__(self, logger_instance = logging.getLogger()):
        self.logger = logger_instance
        self.files = list()
        self.autoname = True
        self.outpath = None
        self.filename = None
        self.modelname = None

    def set_out_path(self, path, filename=None, modelname=None):
        """
        Setter method to set the output path for the IOhandler

        :param path: The requested output path
        :param path_direct: specifies if the path is to be taken
        :return: None
        """
        if modelname is None and filename is None:
            self.autoname = True
            self.outpath = path
        else:
            # Assign name to both model and file name, if only one name is provided
            modelname = filename if modelname is None else modelname
            filename = modelname if filename is None else filename
            # Check if we have string values
            if isinstance(filename, str) and isinstance(modelname, str):
                self.autoname = False
                self.outpath = path
                self.modelname = modelname
                self.filename = filename
            else:
                raise ValueError("Either Modelname or Filename supplied is not a string")


    def load_file(self, path):
        """
        Method to load a sNpfile (Touchstone file) from a given path.
        Loads the file's contents and stores it in an SNpFile class.

        :param path: The path of the file to be loaded
        :return: None
        :raises Exception: if loading the file did not work
        """

        try:
            for actual_path in path:
                ntwk = rf.Network(actual_path)

                #check if file is already loaded -> if so, skip it
                if ntwk.name in [file.name for file in self.files]:
                    self.logger.warning("Warning; file: \"" + ntwk.name + "\" already present, did not load!")
                    continue

                self.logger.info("Opened file: \"" + ntwk.name+"\"")
                self.files.append(ntwk)
        except Exception as e:
            raise e

    def generate_Netlist_2_port(self, parameters, fit_order, fit_type, saturation_table, captype = None):
        """
        Writes an LTSpice Netlist to the path that is stored in the IOhandlers instance variable.

        Will output an LTSpice Netlist for current dependent inductor/capacitor with
        **constant higher order resonances**. That is the higher order resonant circuits **will not** be current/voltage
        dependent in this form of output.


        :param parameters: The Parameters for the model. A Parameters() object containing the model parameters for
            reference file
        :param fit_order: The order of the model i.e. the number of circuits
        :param fit_type: Whether the element is a coil or capacitor
        :param saturation_table: The saturation table for the elements. A dict type object with a key equal to that of
            the parameter in question containing the saturation table as a string
        :param captype: The type of capacitor. Can be GENERIC or MLCC
        :return: None
        """

        out = parameters
        order = fit_order

        match fit_type:
            case constants.El.INDUCTOR:

                if self.autoname:
                    model_name = "L_1"
                else:
                    model_name = self.modelname

                # main element parameters
                L = out['L'].value*config.INDUNIT
                C = out['C'].value*config.CAPUNIT
                R_s = out['R_s'].value
                R_p = out['R_Fe'].value

                lib = '* Netlist for Inductor Model {name} (L={value}H)\n' \
                      '* Including {number} Serially Chained Parallel Resonant Circuits\n*\n'.format(name=model_name,
                                                                                                     value=str(L),
                                                                                                     number=order)
                lib += '.SUBCKT {name} PORT1 PORT2'.format(name=model_name) + '\n*\n'


                ############### HIGHER ORDER ELEMENTS ##################################################################
                for circuit in range(1, order + 1):
                    Cx = out['C%s' % circuit].value*config.CAPUNIT
                    Lx = out['L%s' % circuit].value*config.INDUNIT
                    Rx = out['R%s' % circuit].value
                    node2 = circuit + 1 if circuit < order else 'PORT2'
                    lib += 'C{no} {node1} {node2} '.format(no=circuit, node1=circuit, node2=node2) + str(Cx) + "\n"
                    lib += 'L{no} {node1} {node2} '.format(no=circuit, node1=circuit, node2=node2) + str(Lx) + "\n"
                    lib += 'R{no} {node1} {node2} '.format(no=circuit, node1=circuit, node2=node2) + str(Rx) + "\n"

                ############### MAIN ELEMENT ###########################################################################

                main_res_terminal_port = '1' if order > 0 else 'PORT2'
                lib += 'R_s PORT1 B1 ' + str(R_s) + "\n"
                lib += 'R_p B1 ' + main_res_terminal_port + ' R = limit({lo}, {hi}, {R_Fe} * V(K_Fe))'.format(
                    lo=R_p * 1e-8, hi=R_p * 1e8, R_Fe=R_p) + "\n"

                # B source mimicking the parasitic capacitance
                lib += 'BC PORT1 ' + main_res_terminal_port + ' ' + 'I=-I(BCT)*V(K_C)' + "\n"

                # B source mimicking the main inductor
                lib += 'BL B1 ' + main_res_terminal_port + ' V=V(K_L)*V(LT)' + "\n"

                # 'Test' inductor
                lib += 'L LT 0 ' + str(L) + "\n"
                lib += 'F1 0 LT BL 1' + "\n"

                # 'Test' capacitor
                lib += 'C CT 0 ' + str(C) + '\n'
                lib += 'BCT CT 0 V=V(PORT1)-V(' + main_res_terminal_port + ')' + '\n'

                ############### PROPORTIONALITY TABLES FOR MAIN ELEMENT ################################################
                lib += '* The values for the Current-Inductance-Table can be edited here:' + "\n"
                # proportionality factor for L
                lib += '* current dependent proportionality factor for L' + "\n"
                lib += 'BKL K_L 0 V=table(abs(I(BL)),{table})'.format(table=saturation_table['L']) + "\n"
                # Proportionality factor for C
                lib += '* current dependent proportionality factor for C' + "\n"
                lib += 'BKC K_C 0 V=table(abs(I(BL)),{table})'.format(table=saturation_table['C']) + "\n"
                # proportionality factor for R_Fe
                lib += '* current dependent proportionality factor for R_Fe' + "\n"
                lib += 'BKR K_FE 0 V=table(abs(I(BL)),{table})'.format(table=saturation_table['R_Fe']) + "\n"

                lib += '.ENDS {inductor}'.format(inductor=model_name) + "\n"



            case constants.El.CAPACITOR:

                if self.autoname:
                    model_name = "C_1"
                else:
                    model_name = self.modelname

                # main element parameters
                C = out['C'].value*config.CAPUNIT
                Ls = out['L'].value*config.INDUNIT
                R_s = out['R_s'].value
                R_iso = out['R_iso'].value

                lib = '* Netlist for Capacitor Model {name} (C={value}F)\n' \
                      '* Including {number} Parallely Chained Serial Resonant Circuits\n*\n'.format(name=model_name,
                                                                                                    value=str(C),
                                                                                                    number=order)
                lib += '.SUBCKT {name} PORT1 PORT2'.format(name=model_name) + '\n*\n'

                ############### ACOUSTIC RESONANCE PARAMETERS FOR MLCCs ################################################

                if captype == constants.captype.MLCC:
                    RA = out['R_A'].value
                    LA = out['L_A'].value*config.INDUNIT
                    CA = out['C_A'].value*config.CAPUNIT
                    # current dependent coil for higher order res:
                    lib += 'BL{no} PORT1 NL{node1} '.format(no='A', node1='A') + 'V=V(VL{no})*V(K_L{no})'.format(no='A') + "\n"
                    lib += 'L{no} VL{no} 0 '.format(no='A') + str(LA) + "\n"
                    lib += 'FL{no} 0 VL{no} BL{no} 1'.format(no='A') + "\n"

                    lib += '* current dependent proportionality factor for L{no}'.format(no='A') + "\n"
                    lib += 'BLK{no} K_L{no} 0 V=table(abs(V(PORT1)-V(PORT2)),{table})'.format(no='A',table=saturation_table['L_A']) + "\n"
                    lib += "\n"

                    # current dependent cap for higher order res:
                    lib += 'BC{no} NL{node1} NC{node1} '.format(no='A', node1='A') + 'I=-I(BCT{no})*V(K_C{no})'.format(no='A') + "\n"
                    lib += 'C{no} VC{no} 0 '.format(no='A') + str(CA) + "\n"
                    lib += 'BCT{no} VC{no} 0 '.format(no='A') + 'V=V(NL{node1})-V(NC{node1})'.format(node1='A') + "\n"

                    lib += '* current dependent proportionality factor for C{no}'.format(no='A') + "\n"
                    lib += 'BCK{no} K_C{no} 0 V=table(abs(V(PORT1)-V(PORT2)),{table})'.format(no='A',table=saturation_table['C_A']) + "\n"

                    # current dependent resistor
                    lib += 'R_{no} NC{node1} PORT2 R = limit({lo}, {hi}, {R_x} * V(K_R{no}))'.format(no='A',
                                                                                                     node1='A',
                                                                                                     lo=RA * 1e-8,
                                                                                                     hi=RA * 1e8,
                                                                                                     R_x=RA) + "\n"

                    lib += '* current dependent proportionality factor for R{no}'.format(no='A') + "\n"
                    lib += 'BRK{no} K_R{no} 0 V=table(abs(V(PORT1)-V(PORT2)),{table})'.format(no='A',table=saturation_table['R_A']) + "\n"
                    lib += "\n"

                ############### HIGHER ORDER ELEMENTS ##################################################################
                for circuit in range(1, order + 1):
                    Cx = out['C%s' % circuit].value*config.CAPUNIT
                    Lx = out['L%s' % circuit].value*config.INDUNIT
                    Rx = out['R%s' % circuit].value

                    lib += 'R{no} PORT1 NR{node2} '.format(no=circuit, node2=circuit) + str(Rx) + "\n"
                    lib += 'L{no} NR{node1} NL{node2} '.format(no=circuit, node1=circuit, node2=circuit) + str(
                        Lx) + "\n"
                    lib += 'C{no} NL{node1} PORT2 '.format(no=circuit, node1=order + circuit) + str(Cx) + "\n"

                ############### MAIN ELEMENT ###########################################################################

                # Series resistance
                lib += 'R_s PORT1 LsRs R = limit({lo}, {hi}, {R_s} * V(K_Rs))'.format(lo=R_s * 1e-8, hi=R_s * 1e8,
                                                                                 R_s=R_s) + "\n"
                # Parasitic inductance
                lib += 'L_s LsRs Vcap ' + str(Ls) + "\n"
                # B source mimicking the capacitor
                lib += 'B1 PORT2 Vcap I=I(E1)*V(K_C) ' + "\n"
                # Isolation Resistance
                lib += 'R_iso Vcap PORT2 ' + str(R_iso) + "\n"

                # Test cap (E source taking the voltage over the main capacitance)
                lib += 'E1 VC 0 Vcap PORT2 1 ' + "\n"
                lib += 'C VC 0 ' + str(C) + "\n"

                ############### PROPORTIONALITY TABLES FOR MAIN ELEMENT ################################################
                # Proportionality factor for main element
                lib += '* The values for the Voltage-Capacitance-Table can be edited here:' + "\n"
                lib += 'B2 K_C 0 V=table(abs(V(PORT1)-V(PORT2)),{table}) '.format(table=saturation_table['C']) + "\n"
                # Proportionality factor for series resistance
                lib += '* The values for the Voltage-Resistance-Table can be edited here:' + "\n"
                lib += 'B3 K_Rs 0 V=table(abs(V(PORT1)-V(PORT2)),{table}) '.format(table=saturation_table['R_s']) + "\n"
                # Model closing
                lib += '.ENDS {name}'.format(name=model_name) + "\n"


        ############### OUTPUT #########################################################################################

        if self.autoname:
            # get output folder and path
            out_path = os.path.split(self.outpath)[0]
            dir_name = os.path.normpath(self.outpath).split(os.sep)[-2]
            out_folder = os.path.join(out_path, "fit_results_%s" % dir_name)

            #create the folder; should not be necessary to handle an exception; however folder could be write protected
            try:
                os.makedirs(out_folder, exist_ok = True)
            except Exception:
                raise

            # write LTSpice .lib file
            file = open(os.path.join(out_folder, "LT_Spice_Model_" + dir_name + ".lib"), "w+")
            file.write(lib)
            file.close()
        else:
            out_folder = self.outpath
            # write LTSpice .lib file
            file = open(os.path.join(out_folder, self.filename + ".lib"), "w+")
            file.write(lib)
            file.close()

    def generate_Netlist_2_port_full_fit(self, parameters, fit_order, fit_type, saturation_table, captype=None):
        """
        Writes an LTSpice Netlist to the path that is stored in the IOhandlers instance variable.

        This method **does output fully parametric models**, i.e. the higher order resonant circuits **will be**
        current/voltage dependent as well as the main element.

        :param parameters: The Parameters for the model. A Parameters() object containing the parameters of the
            reference file
        :param fit_order: The order of the model i.e. the number of circuits
        :param fit_type: Whether the element is a coil or capacitor
        :param saturation_table: The saturation table for the elements. A dict type object with a key equal to that of
            the parameter in question containing the saturation table as a string
        :param captype: The type of capacitor. Can be GENERIC or MLCC
        :return: None
        """

        out = parameters
        order = fit_order

        match fit_type:
            case constants.El.INDUCTOR:

                if self.autoname:
                    model_name = "L_1"
                else:
                    model_name = self.modelname

                # parameters for the main elements
                L = out['L'].value*config.INDUNIT
                C = out['C'].value*config.CAPUNIT
                R_s = out['R_s'].value
                R_p = out['R_Fe'].value

                lib = '* Netlist for Inductor Model {name} (L={value}H)\n' \
                      '* Including {number} Serially Chained Parallel Resonant Circuits\n*\n'.format(name=model_name,
                                                                                                     value=str(L), number=order)
                lib += '.SUBCKT {name} PORT1 PORT2'.format(name=model_name) + '\n*\n'

                ############### HIGHER ORDER ELEMENTS ##################################################################
                for circuit in range(1, order + 1):
                    Cx = out['C%s' % circuit].value*config.CAPUNIT
                    Lx = out['L%s' % circuit].value*config.INDUNIT
                    Rx = out['R%s' % circuit].value
                    node2 = circuit + 1 if circuit < order else 'PORT2'

                    #current dependent coil for higher order res:
                    lib += 'BL{no} {node1} {node2} '.format(no=circuit, node1=circuit, node2=node2) + 'V=V(VL{no})*V(K_L{no})'.format(no=circuit) + "\n"
                    lib += 'L{no} VL{no} 0 '.format(no=circuit) + str(Lx) + "\n"
                    lib += 'FL{no} 0 VL{no} BL{no} 1'.format(no=circuit) + "\n"
                    #"test" inductor
                    lib += '* current dependent proportionality factor for L{no}'.format(no=circuit) + "\n"
                    lib += 'BLK{no} K_L{no} 0 V=table(abs(I(BL)),{table})'.format(no=circuit, table=saturation_table['L%s' % circuit]) + "\n"
                    lib += "\n"

                    # current dependent cap for higher order res:
                    lib += 'BC{no} {node1} {node2} '.format(no=circuit, node1=circuit,node2=node2) + 'I=-I(BCT{no})*V(K_C{no})'.format(no=circuit) + "\n"
                    lib += 'C{no} VC{no} 0 '.format(no=circuit) + str(Cx) + "\n"
                    lib += 'BCT{no} VC{no} 0 '.format(no=circuit) + 'V=V({node1})-V({node2})'.format(node1=circuit, node2=node2) + "\n"
                    #"test" cap
                    lib += '* current dependent proportionality factor for C{no}'.format(no=circuit) + "\n"
                    lib += 'BCK{no} K_C{no} 0 V=table(abs(I(BL)),{table})'.format(no=circuit, table=saturation_table['C%s' % circuit]) + "\n"
                    lib += "\n"

                    #current dependent resistor
                    lib += 'R_{no} {node1} {node2} R = limit({lo}, {hi}, {R_x} * V(K_R{no}))'.format(no=circuit, node1=circuit,node2=node2, lo=Rx * 1e-12, hi=Rx * 1e8,R_x=Rx) + "\n"

                    lib += '* current dependent proportionality factor for R{no}'.format(no=circuit) + "\n"
                    lib += 'BRK{no} K_R{no} 0 V=table(abs(I(BL)),{table})'.format(no=circuit, table=saturation_table['R%s' % circuit]) + "\n"
                    lib += "\n"

                ############### MAIN ELEMENT ###########################################################################
                main_res_terminal_port = '1' if order > 0 else 'PORT2'
                lib += 'R_s PORT1 B1 ' + str(R_s) + "\n"
                lib += 'R_p B1 '+main_res_terminal_port+' R = limit({lo}, {hi}, {R_Fe} * V(K_Fe))'.format(lo = R_p * 1e-8, hi = R_p * 1e8, R_Fe = R_p) + "\n"

                # B source mimicking the parasitic capacitance
                lib += 'BC PORT1 ' +main_res_terminal_port+' ' + 'I=-I(BCT)*V(K_C)'+ "\n"

                # B source mimicking the main inductor
                lib += 'BL B1 ' + main_res_terminal_port + ' V=V(K_L)*V(LT)' + "\n"

                # 'Test' inductor
                lib += 'L LT 0 ' + str(L) + "\n"
                lib += 'F1 0 LT BL 1' + "\n"

                # 'Test' capacitor
                lib += 'C CT 0 ' + str(C) + '\n'
                lib += 'BCT CT 0 V=V(PORT1)-V(' + main_res_terminal_port + ')' + '\n'


                ############### PROPORTIONALITY TABLES FOR MAIN ELEMENT ################################################
                lib += '* The values for the Current-Inductance-Table can be edited here:' + "\n"
                #proportionality factor for L
                lib += '* current dependent proportionality factor for L' + "\n"
                lib += 'BKL K_L 0 V=table(abs(I(BL)),{table})'.format(table=saturation_table['L']) + "\n"
                # Proportionality factor for C
                lib += '* current dependent proportionality factor for C' + "\n"
                lib += 'BKC K_C 0 V=table(abs(I(BL)),{table})'.format(table=saturation_table['C']) + "\n"
                #proportionality factor for R_Fe
                lib += '* current dependent proportionality factor for R_Fe' + "\n"
                lib += 'BKR K_FE 0 V=table(abs(I(BL)),{table})'.format(table=saturation_table['R_Fe']) + "\n"

                lib += '.ENDS {inductor}'.format(inductor=model_name) + "\n"


            case constants.El.CAPACITOR:

                if self.autoname:
                    model_name = "C_1"
                else:
                    model_name = self.modelname

                # main element parameters
                C = out['C'].value*config.CAPUNIT
                Ls = out['L'].value*config.INDUNIT
                R_s = out['R_s'].value
                R_iso = out['R_iso'].value

                lib = '* Netlist for Capacitor Model {name} (C={value}F)\n' \
                      '* Including {number} Parallely Chained Serial Resonant Circuits\n*\n'.format(name=model_name,
                                                                                                    value=str(C),
                                                                                                    number=order)
                lib += '.SUBCKT {name} PORT1 PORT2'.format(name=model_name) + '\n*\n'

                ############### HIGHER ORDER ELEMENTS ##################################################################
                for circuit in range(1, order + 1):
                    Cx = out['C%s' % circuit].value*config.CAPUNIT
                    Lx = out['L%s' % circuit].value*config.INDUNIT
                    Rx = out['R%s' % circuit].value
                    node2 = circuit + 1 if circuit < order else 'PORT2'

                    # Voltage dependent coil for higher order res:
                    # B-source mimicking the inductor
                    lib += 'BL{no} PORT1 NL{node1} '.format(no=circuit, node1=circuit) + 'V=V(VL{no})*V(K_L{no})'.format(no=circuit) + "\n"
                    # Test Inductor
                    lib += 'L{no} VL{no} 0 '.format(no=circuit) + str(Lx) + "\n"
                    lib += 'FL{no} 0 VL{no} BL{no} 1'.format(no=circuit) + "\n"
                    # Proportionality factor
                    lib += '* current dependent proportionality factor for L{no}'.format(no=circuit) + "\n"
                    lib += 'BLK{no} K_L{no} 0 V=table(abs(V(PORT1)-V(PORT2)),{table})'.format(no=circuit, table=saturation_table['L%s' % circuit]) + "\n"
                    lib += "\n"

                    # Voltage dependent cap for higher order res:
                    # B source mimicking the capacitor
                    lib += 'BC{no} NL{node1} NC{node1} '.format(no=circuit, node1=circuit) + 'I=-I(BCT{no})*V(K_C{no})'.format(no=circuit) + "\n"
                    # Test cap
                    lib += 'C{no} VC{no} 0 '.format(no=circuit) + str(Cx) + "\n"
                    lib += 'BCT{no} VC{no} 0 '.format(no=circuit) + 'V=V(NL{node1})-V(NC{node1})'.format(node1=circuit) + "\n"
                    # Proportionality factor
                    lib += '* current dependent proportionality factor for C{no}'.format(no=circuit) + "\n"
                    lib += 'BCK{no} K_C{no} 0 V=table(abs(V(PORT1)-V(PORT2)),{table})'.format(no=circuit, table=saturation_table['C%s' % circuit]) + "\n"
                    lib += "\n"

                    # Voltage dependent resistor for higher order res:
                    lib += 'R_{no} NC{node1} PORT2 R = limit({lo}, {hi}, {R_x} * V(K_R{no}))'.format(no=circuit, node1=circuit, lo=Rx * 1e-8, hi=Rx * 1e8,R_x=Rx) + "\n"
                    # Proportionality factor
                    lib += '* current dependent proportionality factor for R{no}'.format(no=circuit) + "\n"
                    lib += 'BRK{no} K_R{no} 0 V=table(abs(V(PORT1)-V(PORT2)),{table})'.format(no=circuit, table=saturation_table['R%s' % circuit]) + "\n"
                    lib += "\n"

                ############### ACOUSTIC RESONANCE FOR MLCCs ###########################################################
                if captype == constants.captype.MLCC:
                    RA = out['R_A'].value
                    LA = out['L_A'].value*config.INDUNIT
                    CA = out['C_A'].value*config.CAPUNIT

                    # Voltage dependent inductance:
                    # B source mimicking the inductor
                    lib += 'BL{no} PORT1 NL{node1} '.format(no='A', node1='A') + 'V=V(VL{no})*V(K_L{no})'.format(no='A') + "\n"
                    # Test inductor
                    lib += 'L{no} VL{no} 0 '.format(no='A') + str(LA) + "\n"
                    lib += 'FL{no} 0 VL{no} BL{no} 1'.format(no='A') + "\n"
                    # Proportionality factor
                    lib += '* current dependent proportionality factor for L{no}'.format(no='A') + "\n"
                    lib += 'BLK{no} K_L{no} 0 V=table(abs(V(PORT1)-V(PORT2)),{table})'.format(no='A', table=saturation_table['L_A']) + "\n"
                    lib += "\n"

                    # Voltage dependent cap:
                    # B source mimicking the cap
                    lib += 'BC{no} NL{node1} NC{node1} '.format(no='A', node1='A') + 'I=-I(BCT{no})*V(K_C{no})'.format(no='A') + "\n"
                    # Test cap
                    lib += 'C{no} VC{no} 0 '.format(no='A') + str(CA) + "\n"
                    lib += 'BCT{no} VC{no} 0 '.format(no='A') + 'V=V(NL{node1})-V(NC{node1})'.format(node1='A') + "\n"
                    # Proportionality factor
                    lib += '* current dependent proportionality factor for C{no}'.format(no='A') + "\n"
                    lib += 'BCK{no} K_C{no} 0 V=table(abs(V(PORT1)-V(PORT2)),{table})'.format(no='A', table=saturation_table['C_A']) + "\n"

                    # Current dependent resistor
                    lib += 'R_{no} NC{node1} PORT2 R = limit({lo}, {hi}, {R_x} * V(K_R{no}))'.format(no='A', node1='A',lo=RA * 1e-8,hi=RA * 1e8,R_x=RA) + "\n"
                    # Proportionality factor
                    lib += '* current dependent proportionality factor for R{no}'.format(no='A') + "\n"
                    lib += 'BRK{no} K_R{no} 0 V=table(abs(V(PORT1)-V(PORT2)),{table})'.format(no='A',table=saturation_table['R_A']) + "\n"
                    lib += "\n"

                ############### MAIN ELEMENT ###########################################################################

                # Series resistance
                lib += 'R_s PORT1 LsRs R = limit({lo}, {hi}, {R_s} * V(K_Rs))'.format(lo=R_s * 1e-8, hi=R_s * 1e8,
                                                                                 R_s=R_s) + "\n"
                # Parasitic inductance
                lib += 'L_s LsRs Vcap ' + str(Ls) + "\n"
                # B source mimicking the capacitor
                lib += 'B1 PORT2 Vcap I=I(E1)*V(K_C) ' + "\n"
                # Isolation Resistance
                lib += 'R_iso Vcap PORT2 ' + str(R_iso) + "\n"

                # Test cap (E source taking the voltage over the main capacitance)
                lib += 'E1 VC 0 Vcap PORT2 1 ' + "\n"
                lib += 'C VC 0 ' + str(C) + "\n"

                ############### PROPORTIONALITY TABLES FOR MAIN ELEMENT ################################################
                # Proportionality factor for main element
                lib += '* The values for the Voltage-Capacitance-Table can be edited here:' + "\n"
                lib += 'B2 K_C 0 V=table(abs(V(PORT1)-V(PORT2)),{table}) '.format(table=saturation_table['C']) + "\n"
                # Proportionality factor for series resistance
                lib += '* The values for the Voltage-Resistance-Table can be edited here:' + "\n"
                lib += 'B3 K_Rs 0 V=table(abs(V(PORT1)-V(PORT2)),{table}) '.format(table=saturation_table['R_s']) + "\n"
                # Model closing
                lib += '.ENDS {name}'.format(name=model_name) + "\n"


        if self.autoname:
            # get output folder and path
            out_path = os.path.split(self.outpath)[0]
            dir_name = os.path.normpath(self.outpath).split(os.sep)[-2]
            out_folder = os.path.join(out_path, "fit_results_%s" % dir_name)

            #create the folder; should not be necessary to handle an exception; however folder could be write protected
            try:
                os.makedirs(out_folder, exist_ok = True)
            except Exception:
                raise

            # write LTSpice .lib file
            file = open(os.path.join(out_folder, "LT_Spice_Model_" + dir_name + ".lib"), "w+")
            file.write(lib)
            file.close()
        else:
            out_folder = self.outpath
            # write LTSpice .lib file
            file = open(os.path.join(out_folder, self.filename + ".lib"), "w+")
            file.write(lib)
            file.close()



    def generate_Netlist_4_port_single_point(self, parametersDM, parametersCM, fit_orderDM, fit_orderCM):

        # define the name of the model here:
        model_name = "CMC_1"

        node1 = "A1"
        node2 = "A2"
        node3 = "B1"
        node4 = "B2"

        lib = '* Netlist for Common Mode Choke Model {name}\n' \

        lib += '.SUBCKT {name} {n1} {n2} {n3} {n4}'.format(name=model_name, n1=node1,n2=node2,n3=node3,n4=node4) + '\n*\n'

        nextnodeA = node1
        nextnodeB = node3

        DMCMnodeA = "DMCMA"
        DMCMnodeB = "DMCMB"


        # DM main resonance
        L = parametersDM['L'].value * config.INDUNIT / 4
        C = parametersDM['C'].value * config.CAPUNIT * 2
        R_s = parametersDM['R_s'].value / 2
        R_p = parametersDM['R_Fe'].value / 2

        lib += 'R_sADM {port1} BDM1A {Rs}'.format(port1=nextnodeA, Rs=R_s) + "\n"
        lib += 'R_pADM BDM1A {MRterminal} {Rp}'.format(MRterminal=DMCMnodeA, Rp=R_p) + "\n"
        lib += 'C_ADM {port1} {MRterminal} {C}'.format(port1=nextnodeA, MRterminal=DMCMnodeA, C=C) + "\n"
        lib += 'L_ADM BDM1A {MRterminal} {L}'.format(MRterminal=DMCMnodeA, L=L) + "\n"
        lib += '\n'

        lib += 'R_sBDM {port1} BDM1B {Rs}'.format(port1=nextnodeB, Rs=R_s) + "\n"
        lib += 'R_pBDM BDM1B {MRterminal} {Rp}'.format(MRterminal=DMCMnodeB, Rp=R_p) + "\n"
        lib += 'C_BDM {port1} {MRterminal} {C}'.format(port1=nextnodeB, MRterminal=DMCMnodeB, C=C) + "\n"
        lib += 'L_BDM BDM1B {MRterminal} {L}'.format(MRterminal=DMCMnodeB, L=L) + "\n"
        lib += '\n'

        lib += 'K{id} {L1} {L2} {K}'.format(id="KDM_main", L1="L_ADM", L2="L_BDM", K="-1") + "\n"

        lib += '\n\n'

        nextnodeA = "MRA" if (fit_orderCM > 0 or fit_orderDM > 0) else node2
        nextnodeB = "MRB" if (fit_orderCM > 0 or fit_orderDM > 0) else node4

        # CM main resonance
        L = parametersCM['L'].value * config.INDUNIT
        C = parametersCM['C'].value * config.CAPUNIT / 2
        R_s = parametersCM['R_s'].value * 2
        R_p = parametersCM['R_Fe'].value * 2

        lib += 'R_sACM {port1} BCM1A {Rs}'.format(port1=DMCMnodeA, Rs=R_s) + "\n"
        lib += 'R_pACM BCM1A {MRterminal} {Rp}'.format(MRterminal=nextnodeA, Rp=R_p) + "\n"
        lib += 'C_ACM {port1} {MRterminal} {C}'.format(port1=DMCMnodeA, MRterminal=nextnodeA, C=C) + "\n"
        lib += 'L_ACM BCM1A {MRterminal} {L}'.format(MRterminal=nextnodeA, L=L) + "\n"
        lib += '\n'

        lib += 'R_sBCM {port1} BCM1B {Rs}'.format(port1=DMCMnodeB, Rs=R_s) + "\n"
        lib += 'R_pBCM BCM1B {MRterminal} {Rp}'.format(MRterminal=nextnodeB, Rp=R_p) + "\n"
        lib += 'C_BCM {port1} {MRterminal} {C}'.format(port1=DMCMnodeB, MRterminal=nextnodeB, C=C) + "\n"
        lib += 'L_BCM BCM1B {MRterminal} {L}'.format(MRterminal=nextnodeB, L=L) + "\n"
        lib += '\n'

        lib += 'K{id} {L1} {L2} {K}'.format(id="KCM_main", L1="L_ACM", L2="L_BCM", K="1") + "\n"

        lib += '\n\n'



        ############### HIGHER ORDER ELEMENTS DM ##################################################################
        for circuit in range(1, fit_orderDM + 1):
            ID = "DM" + str(circuit)

            Cx = (parametersDM['C%s' % circuit].value*2) * config.CAPUNIT
            Lx = (parametersDM['L%s' % circuit].value/4) * config.INDUNIT
            Rx = (parametersDM['R%s' % circuit].value/2)

            n2A = "DM_A_" + str(circuit) if not(circuit == fit_orderDM and fit_orderCM == 0) else node2
            n2B = "DM_B_" + str(circuit) if not(circuit == fit_orderDM and fit_orderCM == 0) else node4

            n1A = nextnodeA
            n1B = nextnodeB


            lib += 'C{no} {node1} {node2} '.format(no=(ID+"A"), node1=n1A, node2=n2A) + str(Cx) + "\n"
            lib += 'L{no} {node1} {node2} '.format(no=(ID+"A"), node1=n1A, node2=n2A) + str(Lx) + "\n"
            lib += 'R{no} {node1} {node2} '.format(no=(ID+"A"), node1=n1A, node2=n2A) + str(Rx) + "\n"

            lib += 'K{id} {L1} {L2} {K}'.format(id=ID, L1="L"+ID+"A", L2="L"+ID+"B",K="-1") + "\n"

            lib += 'C{no} {node1} {node2} '.format(no=(ID+"B"), node1=n1B, node2=n2B) + str(Cx) + "\n"
            lib += 'L{no} {node1} {node2} '.format(no=(ID+"B"), node1=n1B, node2=n2B) + str(Lx) + "\n"
            lib += 'R{no} {node1} {node2} '.format(no=(ID+"B"), node1=n1B, node2=n2B) + str(Rx) + "\n"
            lib += '\n'

            nextnodeA = n2A
            nextnodeB = n2B


        ############### HIGHER ORDER ELEMENTS CM ##################################################################
        for circuit in range(1, fit_orderCM + 1):
            ID = "CM" + str(circuit)

            Cx = (parametersCM['C%s' % circuit].value / 2) * config.CAPUNIT
            Lx = (parametersCM['L%s' % circuit].value) * config.INDUNIT
            Rx = (parametersCM['R%s' % circuit].value * 2)

            n2A = "CM_A_" + str(circuit) if not(circuit == fit_orderCM) else node2
            n2B = "CM_B_" + str(circuit) if not(circuit == fit_orderCM) else node4

            n1A = nextnodeA
            n1B = nextnodeB


            # if circuit < fit_orderCM:
            #     n2A = "CM_A_" + str(circuit + 1)
            #     n2B = "CM_B_" + str(circuit + 1)
            # elif circuit < fit_orderCM and fit_orderDM == 0:
            #     n2A = node2
            #     n2B = node4
            # else:
            #     n2A = DMCMconnectA
            #     n2B = DMCMconnectB



            lib += 'C{no} {node1} {node2} '.format(no=(ID + "A"), node1=n1A, node2=n2A) + str(Cx) + "\n"
            lib += 'L{no} {node1} {node2} '.format(no=(ID + "A"), node1=n1A, node2=n2A) + str(Lx) + "\n"
            lib += 'R{no} {node1} {node2} '.format(no=(ID + "A"), node1=n1A, node2=n2A) + str(Rx) + "\n"

            lib += 'K{id} {L1} {L2} {K}'.format(id=ID, L1="L" + ID + "A", L2="L" + ID + "B", K="1") + "\n"

            lib += 'C{no} {node1} {node2} '.format(no=(ID + "B"), node1=n1B, node2=n2B) + str(Cx) + "\n"
            lib += 'L{no} {node1} {node2} '.format(no=(ID + "B"), node1=n1B, node2=n2B) + str(Lx) + "\n"
            lib += 'R{no} {node1} {node2} '.format(no=(ID + "B"), node1=n1B, node2=n2B) + str(Rx) + "\n"

            lib+= '\n'


            nextnodeA = n2A
            nextnodeB = n2B


        ############### MAIN ELEMENT ###########################################################################



        # Model closing
        lib += '.ENDS {inductor}'.format(inductor=model_name) + "\n"


        ############### OUTPUT #########################################################################################
        # get output folder and path
        out_path = os.path.split(self.outpath)[0]
        dir_name = os.path.normpath(self.outpath).split(os.sep)[-2]
        out_folder = os.path.join(out_path, "fit_results_%s" % dir_name)

        # create the folder; should not be necessary to handle an exception; however folder could be write protected
        try:
            os.makedirs(out_folder, exist_ok=True)
        except Exception:
            raise

        # write LTSpice .lib file
        file = open(os.path.join(out_folder, "LT_Spice_Model_" + dir_name + ".lib"), "w+")
        file.write(lib)
        file.close()

    def generate_Netlist_2_port_single_point(self, parameters, fit_order, fit_type, saturation_table='', captype = None):
        """
        Writes an LTSpice Netlist to the path that is stored in the IOhandlers instance variable.

        Will output an LTSpice Netlist for current dependent inductor/capacitor with
        **constant higher order resonances**. That is the higher order resonant circuits **will not** be current/voltage
        dependent in this form of output.


        :param parameters: The Parameters for the model. A Parameters() object containing the model parameters for
            reference file
        :param fit_order: The order of the model i.e. the number of circuits
        :param fit_type: Whether the element is a coil or capacitor
        :param saturation_table: The saturation table for the elements. A dict type object with a key equal to that of
            the parameter in question containing the saturation table as a string
        :param captype: The type of capacitor. Can be GENERIC or MLCC
        :return: None
        """

        out = parameters
        order = fit_order

        match fit_type:
            case constants.El.INDUCTOR:

                if self.autoname:
                    model_name = "L_1"
                else:
                    model_name = self.modelname

                # main element parameters
                L = out['L'].value*config.INDUNIT
                C = out['C'].value*config.CAPUNIT
                R_s = out['R_s'].value
                R_p = out['R_Fe'].value

                lib = '* Netlist for Inductor Model {name} (L={value}H)\n' \
                      '* Including {number} Serially Chained Parallel Resonant Circuits\n*\n'.format(name=model_name,
                                                                                                     value=str(L),
                                                                                                     number=order)
                lib += '.SUBCKT {name} PORT1 PORT2'.format(name=model_name) + '\n*\n'


                ############### HIGHER ORDER ELEMENTS ##################################################################
                for circuit in range(1, order + 1):
                    Cx = out['C%s' % circuit].value*config.CAPUNIT
                    Lx = out['L%s' % circuit].value*config.INDUNIT
                    Rx = out['R%s' % circuit].value
                    node2 = circuit + 1 if circuit < order else 'PORT2'
                    lib += 'C{no} {node1} {node2} '.format(no=circuit, node1=circuit, node2=node2) + str(Cx) + "\n"
                    lib += 'L{no} {node1} {node2} '.format(no=circuit, node1=circuit, node2=node2) + str(Lx) + "\n"
                    lib += 'R{no} {node1} {node2} '.format(no=circuit, node1=circuit, node2=node2) + str(Rx) + "\n"

                ############### MAIN ELEMENT ###########################################################################

                main_res_terminal_port = '1' if order > 0 else 'PORT2'
                lib += 'R_s PORT1 B1 ' + str(R_s) + "\n"
                lib += 'R_p B1 ' + main_res_terminal_port + ' ' + str(R_p) + "\n"

                # Parasitic capacitance main element
                lib += 'C PORT1 ' + main_res_terminal_port + ' ' + str(C) + "\n"

                # Main inductance
                lib += 'L B1 '+ main_res_terminal_port + ' ' + str(L) + "\n"

                # Model closing
                lib += '.ENDS {inductor}'.format(inductor=model_name) + "\n"



            case constants.El.CAPACITOR:

                if self.autoname:
                    model_name = "C_1"
                else:
                    model_name = self.modelname

                # main element parameters
                C = out['C'].value*config.CAPUNIT
                Ls = out['L'].value*config.INDUNIT
                R_s = out['R_s'].value
                R_iso = out['R_iso'].value

                lib = '* Netlist for Capacitor Model {name} (C={value}F)\n' \
                      '* Including {number} Parallely Chained Serial Resonant Circuits\n*\n'.format(name=model_name,
                                                                                                    value=str(C),
                                                                                                    number=order)
                lib += '.SUBCKT {name} PORT1 PORT2'.format(name=model_name) + '\n*\n'

                ############### MAIN ELEMENT ###########################################################################

                lib += 'R_s PORT1 LsRs ' + str(R_s) + "\n"

                lib += 'L_s LsRs Vcap ' + str(Ls) + "\n"

                lib += 'R_iso Vcap PORT2 ' + str(R_iso) + "\n"

                lib += 'C Vcap PORT2 ' + str(C) + "\n"

                ############### ACOUSTIC RESONANCE PARAMETERS FOR MLCCs ################################################

                ############### HIGHER ORDER ELEMENTS ##################################################################
                for circuit in range(1, order + 1):
                    Cx = out['C%s' % circuit].value*config.CAPUNIT
                    Lx = out['L%s' % circuit].value*config.INDUNIT
                    Rx = out['R%s' % circuit].value

                    lib += 'R{no} PORT1 {node2} '.format(no=circuit, node2=circuit) + str(Rx) + "\n"
                    lib += 'L{no} {node1} {node2} '.format(no=circuit, node1=circuit, node2=order + circuit) + str(
                        Lx) + "\n"
                    lib += 'C{no} {node1} PORT2 '.format(no=circuit, node1=order + circuit) + str(Cx) + "\n"

                # Add acoustic resonance if present
                if captype == constants.captype.MLCC:
                    RA = out['R_A'].value
                    LA = out['L_A'].value*config.INDUNIT
                    CA = out['C_A'].value*config.CAPUNIT

                    lib += 'R{no} PORT1 {node2} '.format(no="A", node2="nA1") + str(RA) + "\n"
                    lib += 'L{no} {node1} {node2} '.format(no="A", node1="nA1", node2="nA2") + str(
                        LA) + "\n"
                    lib += 'C{no} {node1} PORT2 '.format(no="A", node1="nA2") + str(CA) + "\n"


                lib += '.ENDS {name}'.format(name=model_name) + "\n"


        ############### OUTPUT #########################################################################################
        if self.autoname:
            # get output folder and path
            out_path = os.path.split(self.outpath)[0]
            dir_name = os.path.normpath(self.outpath).split(os.sep)[-2]
            out_folder = os.path.join(out_path, "fit_results_%s" % dir_name)

            #create the folder; should not be necessary to handle an exception; however folder could be write protected
            try:
                os.makedirs(out_folder, exist_ok = True)
            except Exception:
                raise

            # write LTSpice .lib file
            file = open(os.path.join(out_folder, "LT_Spice_Model_" + dir_name + ".lib"), "w+")
            file.write(lib)
            file.close()
        else:
            out_folder = self.outpath
            # write LTSpice .lib file
            file = open(os.path.join(out_folder, self.filename + ".lib"), "w+")
            file.write(lib)
            file.close()

    def export_parameters(self, param_array, order, fit_type, captype = None):
        """
        Method to output the obtained model parameters as an .xlsx file to the directory in IOhandlers output path.

        :param param_array: An array containing Parameters() type objects for each file
        :param order: The order of the model, i.e. the number of resonance circuits
        :param fit_type: Whether the model is for a coil or capacitor
        :param captype: The type of capacitor. Can be GENERIC or MLCC
        :return: None
        """

        out_dict = {}

        #write the main resonance parameters to the dict
        match fit_type:
            case constants.El.INDUCTOR:
                R_s_list = []
                R_Fe_list =[]
                L_list = []
                C_list = []
                for param_set in param_array:
                    R_s_list.append(param_set['R_s'].value)
                    R_Fe_list.append(param_set['R_Fe'].value)
                    L_list.append(param_set['L'].value*config.INDUNIT)
                    C_list.append(param_set['C'].value*config.CAPUNIT)

                out_dict['R_s'] = R_s_list
                out_dict['R_Fe'] = R_Fe_list
                out_dict['L'] = L_list
                out_dict['C'] = C_list


            case constants.El.CAPACITOR:
                R_s_list = []
                R_Iso_list = []
                L_list = []
                C_list = []
                for param_set in param_array:
                    R_s_list.append(param_set['R_s'].value)
                    R_Iso_list.append(param_set['R_iso'].value)
                    L_list.append(param_set['L'].value*config.INDUNIT)
                    C_list.append(param_set['C'].value*config.CAPUNIT)

                out_dict['R_s'] = R_s_list
                out_dict['R_iso'] = R_Iso_list
                out_dict['L'] = L_list
                out_dict['C'] = C_list

                if captype == constants.captype.MLCC:
                    R_A_list = []
                    L_A_list = []
                    C_A_list = []
                    for param_set in param_array:
                        R_A_list.append(param_set['R_A'].value)
                        L_A_list.append(param_set['L_A'].value*config.INDUNIT)
                        C_A_list.append(param_set['C_A'].value*config.CAPUNIT)
                    out_dict['R_A'] = R_A_list
                    out_dict['L_A'] = L_A_list
                    out_dict['C_A'] = C_A_list





        for key  in range(1,order+1):
            #generate key numbers and empty lists for the parameters
            C_key = "C%s" % key
            L_key = "L%s" % key
            R_key = "R%s" % key
            w_key = "w%s" % key
            BW_key = "BW%s" % key

            clist = []
            llist = []
            rlist = []
            wlist = []
            bwlist = []

            #iterate through parameter sets
            for param_set in param_array:
                clist.append(param_set[C_key].value*config.CAPUNIT)
                llist.append(param_set[L_key].value*config.INDUNIT)
                rlist.append(param_set[R_key].value)
                wlist.append(param_set[w_key].value*config.FUNIT)
                bwlist.append(param_set[BW_key].value*config.FUNIT)

            out_dict[C_key] = clist
            out_dict[L_key] = llist
            out_dict[R_key] = rlist
            out_dict[w_key] = wlist
            out_dict[BW_key] = bwlist

        #write parameters to a pandas dataframe and transpose
        data_out = pd.DataFrame(out_dict)
        # data_out.transpose()
        if self.autoname:
            out_path = os.path.split(self.outpath)[0]
            dir_name = os.path.normpath(self.outpath).split(os.sep)[-2]
            out_folder = os.path.join(out_path, "fit_results_%s" % dir_name)
            try:
                os.makedirs(out_folder, exist_ok = True)
            except Exception:
                raise

            data_out.to_excel(os.path.join(out_folder, "Parameters_" + dir_name + ".xlsx"))

        else:
            out_folder = self.outpath
            data_out.to_excel(os.path.join(out_folder, "Parameters_" + self.filename + ".xlsx"))

    def output_plot(self, freq, z21, mag, ang, mdl, filename):
        """
        Method to output a Bode-plot and a linear difference plot of the model.

        :param freq: The frequency vector
        :param z21: The measured impedance data
        :param mag: The measured, smoothed magnitude data
        :param ang: The measured, smoothed phase data
        :param mdl: The model data (complex)
        :param filename: The name of the file that the plot will be made for
        :return: None
        """
        if self.autoname:
            out_path = os.path.split(self.outpath)[0]
            dir_name = os.path.normpath(self.outpath).split(os.sep)[-2]
            out_folder = os.path.join(out_path, "fit_results_%s" % dir_name)
            plot_folder = os.path.join(out_folder, "plots")

            try:
                os.makedirs(out_folder, exist_ok=True)
            except Exception:
                raise

            try:
                os.makedirs(plot_folder, exist_ok=True)
            except Exception:
                raise
        else:
            plot_folder = self.outpath


        title = filename
        # fig = plt.figure(figsize=(20, 20))
        fig, ax = plt.subplots(nrows=2,ncols=1)
        fig.set_figheight(20)
        fig.set_figwidth(20)
        #file_title = get_file_path.results + '/03_Parameter-Fitting_' + file_name + "_" + mode
        # plt.subplot(211)
        fig = plt.gcf()
        fig.suptitle(str(title), fontsize=25, fontweight="bold")
        ax[0].set_xscale('log')
        ax[0].set_yscale('log')
        ax[0].set_xlim([min(freq), max(freq)])
        ax[0].set_ylabel('Magnitude in \u03A9', fontsize=16)
        ax[0].set_xlabel('Frequency in Hz', fontsize=16)
        ax[0].grid(True, which="both")
        ax[0].tick_params(labelsize=16)
        ax[0].tick_params(labelsize=16)
        ax[0].plot(freq, abs(z21), 'r', linewidth=3, alpha=0.33, label='Measured Data')
        ax[0].plot(freq, mag, 'r', linewidth=3, alpha=1, label='Filtered Data')
        # Plot magnitude of model in blue
        ax[0].plot(freq, abs(mdl), 'b--', linewidth=3, label='Model')
        ax[0].legend(fontsize=16)
        #Phase
        curve = np.angle(z21, deg=True)
        ax[1].set_xscale('log')
        ax[1].set_xlim([min(freq), max(freq)])
        ax[1].set_ylabel('Phase in °', fontsize=16)
        ax[1].set_xlabel('Frequency in Hz', fontsize=16)
        ax[1].grid(True, which="both")
        ax[1].set_yticks(np.arange(45 * (round(min(curve) / 45)), 45 * (round(max(curve) / 45)) + 1, 45.0))
        ax[1].tick_params(labelsize=16)
        ax[1].tick_params(labelsize=16)
        ax[1].plot(freq, np.angle(z21, deg=True), 'r', linewidth=3, zorder=-2, alpha=0.33, label='Measured Data')
        ax[1].plot(freq, ang, 'r', linewidth=3, zorder=-2, alpha=1, label='Filtered Data')
        #   Plot Phase of model in magenta
        ax[1].plot(freq, np.angle(mdl, deg=True), 'b--', linewidth=3, label='Model', zorder=-1)
        #plt.scatter(resonances_pos, np.zeros_like(resonances_pos) - 90, linewidth=3, color='green', s=200, marker="2",
        #            label='Resonances')
        ax[1].legend(fontsize=16)

        #may be obsolete
        plt.savefig(os.path.join(plot_folder, "Bode_plot_" + filename + ".png"), dpi = 300)

        if constants.SHOW_BODE_PLOTS:
            plt.show()
        else:
            plt.close(fig)

        #Diffplot
        if constants.OUTPUT_DIFFPLOTS:
            diff_data = abs(mdl)-abs(z21)
            diff_data_percent = (diff_data/abs(z21))*100
            title = filename + " (Model-Measurement)/Measurement in %"
            fig = plt.figure(figsize=(20, 20))
            plt.plot(freq, diff_data_percent, 'r', linewidth=3, alpha=1)
            plt.title((title), fontsize=25, fontweight="bold")
            plt.xscale('log')
            plt.yscale('linear')
            plt.xlim([min(freq), max(freq)])
            plt.ylabel('Error in %', fontsize=16)
            plt.xlabel('Frequency in Hz', fontsize=16)
            plt.xticks(fontsize=16)
            plt.yticks(fontsize=16)
            plt.grid(True, which="both")
            # plt.plot(freq, diff_data_percent, 'r', linewidth=3, alpha=1)
            plt.savefig(os.path.join(plot_folder, "Diff_plot_" + filename + ".png"))

        if constants.SHOW_BODE_PLOTS:
            plt.show()
        else:
            plt.close(fig)


