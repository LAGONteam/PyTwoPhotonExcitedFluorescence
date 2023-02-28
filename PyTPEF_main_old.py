"""
Created on 5/1/2023
@author: Jonathan DANIEL, jonathan.daniel@u-bordeaux.fr
"""

import logging
import time

logging.basicConfig(level=logging.INFO,
                    filename="PyTPEF.log",
                    filemode="w",
                    format='%(asctime)s - %(levelname)s - %(message)s')

logging.info("Main import => ok !")


SIMULATION_FOR_DEBUG= False  #Set True only for code test and debug, for real measure set False
                            # do not forget to do the same for configmanagement in Ressources_scripts

import json
import os
import sys

from PyQt6.QtWidgets import QApplication, QWidget, QListWidgetItem, QFileDialog, QMessageBox

pop_up = QApplication(sys.argv)
root = QFileDialog.getExistingDirectory(caption="Please select the folder to save datas")
CUR_DIR = os.path.dirname(__file__)

temp = f"{CUR_DIR}/Ressources/save_directory.json"
root_for_logging = f"{CUR_DIR}/Ressources"

with open(temp, "w") as w:
    json.dump(root, w)

from time import sleep
import numpy as np
import pyqtgraph as pg
from pyqtgraph import PlotWidget, plot
from PyQt6.QtCore import Qt, QTimer, QProcess, QMutex
from PySide6.QtCore import QObject, QThread, Signal
from PyQt6 import uic
from PyQt6.QtGui import QPixmap
from PyQt6.uic.Compiler.qtproxies import QtGui

if SIMULATION_FOR_DEBUG == False:
    from Ressources_scripts.coherent_laser import Chameleon as fs_laser
    from Ressources_scripts import talk_to_spectro
    from Ressources_scripts import talk_to_elliptec_devices as tuner
    import nidaqmx
else:
    from Ressources_scripts.Laser_Simulation import Chameleon as fs_laser
    from Ressources_scripts import Spectro_Simulation as talk_to_spectro
    from Ressources_scripts import Motor_Simulation as tuner
    import Ressources_scripts.Photodiode_Simulation as nidaqmx

from Ressources_scripts.configmanagement import sample_data, save, parrallel_execution, create_json, chrono, get_samples, \
    create_samples, get_path
from Ressources_scripts.tpa_calculation import Read_Data_File as rdf
from Ressources_scripts.configmanagement import Config as cfg
from Ressources_scripts.Data_processing_tpa_calculation import Read_Data_File as rdf_2
from Ressources_scripts.Data_processing import save as sdp
from Ressources_scripts import References_data
from statistics import mean
from pathlib import Path

if SIMULATION_FOR_DEBUG:
    logging.warning("PyTPEF executed in Simulation mode !")

ANGLE_POWER_DICT={}

class Measure(QObject):
    """
    This class is used to do the measurements outside the gui interface to avoid freezing of the application.
    it emits signals at different steps of the measurement and require arguments to start.

    list_of_square_dark_corr_power, list_of_fully_corrected_area, fit, sample, wavelength

    signals:

        signal_to_plot: (array1, array2, float1, int1, int2, float2) send emission spectrum from the spectro
            (array1: wavelength, array2: intensity), P² (float1), F (int1), measure n° (int2), and corresponding
            power from angle_power dictionary (float2) to the plotwidget FluoIntensitySpectrum.

        plot_to_save: (array1, array2, int, str) save the dark emission spectrum (array1: wavelength, array2:intensity)
            to the sample (str) at specific wavelength (int) as a png file.

        measure_finished: (dict1, dict2, dict3, dict4) send sample informations (dict1), experimental data (dict2)
            processed data (dict3) and sample root (dict4).

        finished: end of the measure.

        clear_plot: clear the plotwidgets FluoIntensitySpectrum, FPSpectrum_2 and S2Spectrum.

        save_fluo_vs_power: (int1, list1, list2, str, int2) send number of measure (int1), data from spectro
            (list1: wavelength,list2: intensity), sample name (str) and excitation wavelength (int2) to plot and save
            the emission against power figure as png.


        save_quad: (list1, list2, float, str, int) send P² (list1), F (list2), fit value (float, from log(F)/log(P²))
            for sample (str) at specific excitation wavelength (int) to plot and save the quadraticity check as a png.

        progress_bar_status: (int), send the progression step to gui to help user knowing the application is not freeze.

        signal_to_pop_up: send the pop-up of the end of measure to the user.

        enable_widgets: gives the access of key buttons

        send_data_to_user: (str) send str informations to the user using the ls_terminal_information_for_users_measurement.

    """
    signal_to_plot= Signal(object, object, object, object, object, object)
    plot_to_save= Signal(object, object, object, object)
    measure_finished= Signal(object, object, object, object)
    finished = Signal()
    clear_plot= Signal()
    save_fluo_vs_power= Signal(object, object, object, object, object)
    save_quad= Signal(object, object, object, object, object)
    progress_bar_status=Signal(object)
    signal_to_pop_up= Signal()
    enable_widgets= Signal()
    send_data_to_user = Signal(object)


    def __init__(self, sample, wavelength, integration_time, dictionary_of_sample_info,
                 dictionary_of_sample_processed_data, dictionary_of_sample_experimental_data, fluorescence_ref_dye,
                 dictionary_of_sample_root, ccd_correction_dye, ccd_correction_factor, maximum_power_requested,
                 minimum_power_requested, power_pitch, number_of_scan, new_measure, data_for_new_measure,
                 current_emission_filter, opo_measure, real_opo_wavelength):

        """
        :param sample: str, name of the sample
        :param wavelength: int, excitation wavelength
        :param integration_time: int, integration time send to the spectro for measurement
        :param dictionary_of_sample_info: dict, dictionary containing sample informations (phi, concentration, solvent..)
        :param dictionary_of_sample_processed_data: dict, dictionary containing processed data
        :param dictionary_of_sample_experimental_data: dict, dictionary containing experimental data
        :param fluorescence_ref_dye: str, name of the fluorescent dye used as reference
        :param dictionary_of_sample_root: dict, dictionary containing the sample root
        :param ccd_correction_dye: str, name of the fluorescent dye used to make CCD correction
        :param ccd_correction_factor: array, correction factor for each wavelength of the CCD, size and corresponding
        wavelength depend of the dye used for the CCD correction
        :param maximum_power_requested: int, maximum value of power for measuremnt set by user
        :param minimum_power_requested: int, minimum value of power for measuremnt set by user
        :param power_pitch: int, corresponds to the number of total measure -2
        :param number_of_scan: int, number of accumulation for spectro
        :param new_measure: bool, is the measurement another measure of the sample at this wavelength
        :param data_for_new_measure: dict, if new_measure == True , this dict contains the old data to update
        :param current_emission_filter: str, name of the short pass edge filter (to remove laser signal)
        :param opo_measure: bool, allow measure to know is opo is active
        :param real_opo_wavelength: if opo_measure == True, gives the wavelength of the opo cavity not the pump (800 nm)
        """

        super().__init__()
        """Variables"""
        self.sample= sample
        self.integration_time= integration_time
        self.fluorescence_ref_dye= fluorescence_ref_dye
        self.ccd_correction_dye= ccd_correction_dye
        self.ccd_correction_factor= ccd_correction_factor
        self.maximum_power_requested= maximum_power_requested
        self.minimum_power_requested= minimum_power_requested
        self.power_pitch= power_pitch
        self.number_of_scan= number_of_scan
        self.new_measure= new_measure
        self.sample_info= dictionary_of_sample_info
        self.experimental_data= dictionary_of_sample_experimental_data
        self.processed_data= dictionary_of_sample_processed_data
        self.fluorescence_ref_dye= fluorescence_ref_dye
        self.root_dict = dictionary_of_sample_root
        self.data_for_new_measure=data_for_new_measure
        self.current_emission_filter= current_emission_filter
        self.wavelength=wavelength
        if opo_measure :
            self.real_wavelength=real_opo_wavelength
        else:
            self.real_wavelength=wavelength

        """Spectro"""
        self.spectro= talk_to_spectro

        logging.info("Initialisation of Measure => ok !")

    def run_measure(self):
        """
        Start the measurement of the two-photon excited fluorescence (TPEF)
        """

        self.progress_bar_status.emit(5)
        self.send_data_to_user.emit("Start measuring.")

        logging.info(f"self.sample_info: {self.sample_info.keys()}")

        sample_info, experimental_data, processed_data, root_dict= self.measure_step_1(
                  maximum_power_requested= self.maximum_power_requested,
                  minimum_power_requested= self.minimum_power_requested,
                  power_pitch= self.power_pitch,
                  wavelength= self.real_wavelength,
                  sample= self.sample,
                  integration_time= self.integration_time,
                  number_of_scans= self.number_of_scan,
                  new_measure= self.new_measure
                  )
        save().fill_angle_power_dict_for_measure(data=ANGLE_POWER_DICT)

        logging.info("Fill_angle_power_dict_for_measure saved => ok !")

        self.measure_finished.emit(sample_info, experimental_data, processed_data, root_dict)
        self.signal_to_pop_up.emit()
        self.enable_widgets.emit()
        self.send_data_to_user.emit("Ready for next measure.")
        self.progress_bar_status.emit(100)
        self.finished.emit()

    def get_power_data(self, maximum_power_requested, minimum_power_requested, power_pitch, wavelength):
        """
        Get data from talk_to_elliptec.
        :param maximum_power_requested: int
        :param minimum_power_requested: int
        :param power_pitch: int
        :param wavelength: int
        :return: list (angles), list (power), dict (power_values)
        """

        """Get the corresponding angles for requested max, min & steps power at selected wavelength."""
        angle_min, angle_max, dictionary = tuner.RotationMount().determine_angle_for_power(
            power_max=maximum_power_requested,
            power_min=minimum_power_requested,
            number_of_measure=power_pitch, wavelength=wavelength)

        logging.info("Determine_angle_for_power => ok !")

        """Store angles and corresponding power as lists."""
        list_of_angles, list_of_power, all_power_values = tuner.RotationMount().angles_for_measure(angle_min=angle_min,
                                                                                                   angle_max=angle_max,
                                                                                                   dictionary=dictionary,
                                                                                                   number=power_pitch)

        logging.info("Angles_for_measure => ok !")

        for n in range(len(list_of_angles)):
            angle=list_of_angles[n]
            power=list_of_power[n]
            data=[angle, power]
            if n == 0:
                ANGLE_POWER_DICT.update({f"{wavelength}": {f"{n}":data}})
            else:
                ANGLE_POWER_DICT[f"{wavelength}"].update({f"{n}":data})


        return list_of_angles, list_of_power, all_power_values

    def measure_step_1(self, maximum_power_requested, minimum_power_requested, power_pitch, wavelength, sample,
                       integration_time, number_of_scans, new_measure):
        """
        Check if the measure is a new measure or not and send data to measure_step_2
        :param maximum_power_requested: int
        :param minimum_power_requested: int
        :param power_pitch: int
        :param wavelength: int
        :param sample: str
        :param integration_time: int
        :param number_of_scans: int
        :param new_measure: bool
        :return: dict (sample_info), dict (experimental_data), dict (processed_data), dict (root)
        """
        """Set the laser at desired wavelength."""

        logging.info("Starting measure step 1.")

        fs_laser().setWavelengthBlocking(int(self.wavelength))

        logging.info(f"Laser tuned to {self.wavelength} nm => ok !")

        self.send_data_to_user.emit(f"Tuning the wavelength of the laser to {wavelength} nm.")

        logging.info(f"Max. Power = {maximum_power_requested}.")
        logging.info(f"Min. Power = {minimum_power_requested}.")
        logging.info(f"Pitch power = {power_pitch}.")

        self.send_data_to_user.emit(f"Max. Power = {maximum_power_requested} \n Min. Power = {minimum_power_requested} \n Pitch power = {power_pitch}")
        if new_measure == False:
            """Get list of angles and corresponding power to measure."""
            list_of_angles, list_of_power, all_power_values = self.get_power_data(
                maximum_power_requested=maximum_power_requested, minimum_power_requested=minimum_power_requested,
                power_pitch=power_pitch, wavelength=wavelength)

            logging.info("get_power_data => ok !")

            self.progress_bar_status.emit(15)
        elif new_measure == True:
            list_of_angles=[]
            all_power_values={}
            for i in self.data_for_new_measure:
                temp_angle, temp_power=save().extract_angle_power_from_dict(wavelength=wavelength, i=i)
                list_of_angles.append(temp_angle)
                all_power_values.update({f"{temp_angle}": f"{temp_power}"})

                logging.info(f"New measure of #{i}, angle: {temp_angle}, power: {temp_power}.")

                self.progress_bar_status.emit(15)
        self.measure_step_2(sample=sample, integration_time=integration_time, number_of_scans=number_of_scans,
                            list_of_angles=list_of_angles, wavelength=wavelength, all_power_values=all_power_values,
                            new_measure=new_measure)
        self.progress_bar_status.emit(90)

        logging.info("Measure step 1 ended with success => ok !")

        return self.sample_info, self.experimental_data, self.processed_data, self.root_dict

    def measure_step_2(self, sample, integration_time, number_of_scans, list_of_angles, wavelength, all_power_values,
                       new_measure):
        """
        Perform the dark measurement, send data to measure_step_3 and update dictionary data.
        :param sample: str
        :param integration_time: int
        :param number_of_scans: int
        :param list_of_angles: list
        :param wavelength: int
        :param all_power_values: dict
        :param new_measure: bool
        """
        """Create temporary lists"""

        logging.info("Starting measure step 2.")

        list_of_fully_corrected_area = []
        list_of_square_dark_corr_power = []
        temporary_wavelength_list = []
        temporary_intensity_list = []
        print(sample, "is analyse")

        logging.info(f"Analysing {sample}.")
        logging.info("Analysing the dark.")

        self.send_data_to_user.emit(f"{sample} is analyse.")
        self.send_data_to_user.emit('Measuring the dark.')
        """Dark measuring"""
        talk_to_spectro.Maya(integration_time)

        logging.info(f"Integration time :{integration_time} ms, applyied to Spectro => ok !")

        sleep(0.1)
        x_dark, y_dark, power_dark, raw_fluo = parrallel_execution(integration_time=integration_time).run(
            scan=1, intensity_dark=0)

        logging.info("Parrallel_execution done => ok !")
        logging.info(f"Dark measured. Power dark= {power_dark}.")

        self.send_data_to_user.emit(f"Dark measured with success.")
        self.progress_bar_status.emit(20)
        """Get the mean value of the power"""
        power_dark_mean = mean(power_dark)

        logging.info(f"Mean power for dark= {power_dark_mean}.")

        self.send_data_to_user.emit(f"The mean value of the dark is: {round(float(power_dark_mean),6)}.")
        temporary_wavelength_list.append(x_dark)
        temporary_intensity_list.append(y_dark)

        logging.info("Experimental dark data saved => ok !")

        self.send_data_to_user.emit("Save experimental dark data, ok !")
        """Create and save the dark as a figure"""
        self.plot_to_save.emit(x_dark, y_dark, wavelength, sample)
        self.clear_plot.emit()
        self.send_data_to_user.emit("Dark plot with success.")
        i = 0
        """For each angles defined in the list_of_angles list, measure the fluorescence & power, correct the 
        fluorescence signal, then return the data as lists"""
        """Open the shutter of the laser."""
        self.send_data_to_user.emit("Opening the shutter of the laser !")
        fs_laser().openShutterBlocking()

        logging.info("Laser shutter open => Ok !")

        self.progress_bar_status.emit(25)
        """Wait 15 sec. to let the laser being stable."""
        self.send_data_to_user.emit("Waiting 15s. for laser stabilization.")
        time_shift=0
        for unit_time in range(15):
            sleep(1)
            self.progress_bar_status.emit(26+time_shift)
            if unit_time == 5:
                self.send_data_to_user.emit("Waiting 10s. for laser stabilization.")
            elif unit_time == 10:
                self.send_data_to_user.emit("Waiting 5s. for laser stabilization.")
            time_shift+=1

        progress = 20/int(len(list_of_angles))
        n=1

        for angle in list_of_angles:
            temporary_intensity_list, temporary_wavelength_list, i, list_of_fully_corrected_area, list_of_square_dark_corr_power, full_corrected_area = self.measure_step_3(
                wavelength=wavelength, sample=sample, i=i, angle=angle, all_power_values=all_power_values,
                y_dark=y_dark, power_dark_mean=power_dark_mean,
                temporary_intensity_list=temporary_intensity_list,
                temporary_wavelength_list=temporary_wavelength_list,
                list_of_fully_corrected_area=list_of_fully_corrected_area,
                list_of_square_dark_corr_power=list_of_square_dark_corr_power, new_measure=new_measure,
                integration_time=integration_time, number_of_scans=number_of_scans)
            self.progress_bar_status.emit(40+(progress*n))
            n+=1
        self.progress_bar_status.emit(60)
        self.send_data_to_user.emit("Closing the shutter of the laser.")
        """Close the shutter of the laser."""
        fs_laser().closeShutterBlocking()

        logging.info("Laser shutter close => ok !")
        logging.info("End of measure.")

        self.send_data_to_user.emit("End of measure.")

        logging.info(f"Type of wavelength list for plot: {type(temporary_wavelength_list)}.")
        logging.info(f"Type of intensity list for plot: {type(temporary_intensity_list)}.")
        logging.info(f"Number of measure= {i-1}.")

        self.send_data_to_user.emit("Generating data for plotting.")
        """plot figure_intensity_vs_power"""
        self.save_fluo_vs_power.emit(i, temporary_wavelength_list,temporary_intensity_list, sample, wavelength)

        logging.info("figure fluo_vs_power generated => ok !")

        self.progress_bar_status.emit(65)

        """Calculate the parameters of the fit of F vs P² and return them for saving."""
        area, correlation_coeff, slope, coeff, log_square_fluo, log_square_power, fit = rdf()._Quad_Log(list_of_fully_corrected_area,
                                                                               list_of_square_dark_corr_power)
        logging.info("Quadraticty calculation => ok !")

        self.progress_bar_status.emit(70)
        state = self.sample_info[f"{sample}"]["state"]

        logging.info(f"Sate of {sample}: {state}, {self.sample_info.keys()}.")

        """Save experimental datas"""
        file_exp = save().save_exp(dictionary=self.experimental_data[f"{sample}"], sample=sample, wavelength=wavelength)
        self.send_data_to_user.emit("Data saved.")

        logging.info("Data saved.")

        """Plot and saved the quadraticity as a figure."""
        self.save_quad.emit(list_of_square_dark_corr_power, list_of_fully_corrected_area, fit, sample, wavelength)
        self.progress_bar_status.emit(75)

        logging.info(f"Quadraticity of {sample} at {wavelength} nm saved => ok !")

        """Save processed data."""
        sigma2 = ""
        sigma2_phi = ""
        if self.sample_info[f"{sample}"]["state"] == True:
            self.create_dictionary_processed_data(sample=sample, wavelength=wavelength, slope=slope,
                                                  coeff=coeff, log_square_fluo=log_square_fluo,
                                                  log_square_power=log_square_power, fit=fit,
                                                  full_corrected_area=area,
                                                  sigma2_phi=sigma2_phi, sigma2=sigma2)
            print("////"*100, "State = True")

            logging.info("Creation of dictionary_processed_data => ok !")

        else:
            self.update_dictionary_processed_data(sample=sample, wavelength=wavelength, slope=slope,
                                                  coeff=coeff, log_square_fluo=log_square_fluo,
                                                  log_square_power=log_square_power, fit=fit,
                                                  full_corrected_area=area,
                                                  sigma2_phi=sigma2_phi, sigma2=sigma2)

            logging.info("Update of dictionary_processed_data => ok!")

        self.sample_info[f"{sample}"]["state"] = False
        print("state of dict is : ", self.sample_info[f"{sample}"]["state"])

        logging.info(f"State of self.sample_info[{sample}][state]: False.")

        self.progress_bar_status.emit(80)
        sigma2_phi, sigma2 = rdf()._Calcul_TPA(sample_info=self.sample_info,
                                               processed_data=self.processed_data,
                                               sample_name=sample,
                                               reference_name=self.fluorescence_ref_dye,
                                               wavelength=wavelength)

        logging.info("TPA calculation => ok !")

        self.update_dictionary_processed_data(sample=sample, wavelength=wavelength, slope=slope, coeff=coeff,
                                              log_square_fluo=log_square_fluo,
                                              log_square_power=log_square_power, fit=fit,
                                              full_corrected_area=area, sigma2_phi=sigma2_phi,
                                              sigma2=sigma2)

        logging.info("Update of dictionary_processed_data => ok !")

        #logging.info(f"self.processed_data: {self.processed_data}.")

        sleep(0.1)
        self.progress_bar_status.emit(85)
        file_process = save().save_process(dictionary=self.processed_data[sample], sample=sample, state=False)

        logging.info(f"Processed data of {sample} saved as a json file.")

        save().save_process_data_to_excel(file=file_process, sample=sample)

        logging.info(f"Processed data of {sample} saved as an excel file.")

        self.root_dict.update({sample: file_process})

        logging.info("Measure step 2 ended with success => ok !")
        logging.info(f"Measure of {sample} is finished.")

    def measure_step_3(self, wavelength, sample, i, angle, all_power_values, y_dark, power_dark_mean,
                       temporary_wavelength_list, temporary_intensity_list, list_of_fully_corrected_area,
                       list_of_square_dark_corr_power, new_measure, integration_time, number_of_scans):
        """
        Set the desire power of the laser then measure the emission spectrum & the power of the laser
        :param wavelength: int
        :param sample: str
        :param i: int
        :param angle: int
        :param all_power_values: dict
        :param y_dark: array
        :param power_dark_mean: float
        :param temporary_wavelength_list: list
        :param temporary_intensity_list: list
        :param list_of_fully_corrected_area: list
        :param list_of_square_dark_corr_power: list
        :param new_measure: bool
        :param integration_time: int
        :param number_of_scans: int
        :return: list (intensity), list (wavelength), int (measure n°), list (list of F), list (P²), int (F)
        """

        logging.info("Starting measure step 3.")

        """Set the half-wave plate to the desired angle to get desired power."""
        self.send_data_to_user.emit(f"Turning the half wave plate to {angle}°.")
        tuner.RotationMount().spin_to_position(angle)

        logging.info("Spin_to_position => ok !")

        sleep(0.2)
        true_power = all_power_values[angle]

        logging.info(f"Going to angle= {angle} °.")
        logging.info(f"Real power= {true_power} mW.")
        logging.info(f"Current wavelength= {wavelength} nm.")

        self.send_data_to_user.emit(f"Rotation ok !\nPower is {true_power} mW and wavelength is {wavelength} nm.")
        """Record the corrected fluorescence & power."""
        spectro_wavelength, spectro_intensity, photodiode_power, raw_fluorescence = parrallel_execution(
            integration_time).run(scan=number_of_scans, intensity_dark=y_dark)

        logging.info("Parrallel_execution => ok !")

        """Calculated the mean of power."""
        try:
            photodiode_power_mean = mean(photodiode_power)
        except:
            photodiode_power_mean = 1

            logging.error(f"Can't calculate mean({photodiode_power}), photodiode_power_mean is set to 1.")

        logging.info(f"Measurement of {sample} at angle position {angle} ({true_power} mW) done.")

        self.send_data_to_user.emit(f"Measurement of {sample} at angle position {angle} ({round(float(true_power), 5)} mW) done.\nCalculating and storing dark corrected intensity and power values.")
        """Calculate and store dark corrected intensity and power values"""
        intensity_dark_corr = spectro_intensity

        logging.info(f"The type of intensity dark corr is: {type(intensity_dark_corr)}.")

        self.send_data_to_user.emit("Converting and correcting fluorescence data.")
        """Convert the raw and corrected fluorescence from table to list."""
        spectro_intensity.tolist()
        raw_fluorescence.tolist()

        logging.info(f"The type of data of spectro intensity is: {type(spectro_intensity)}.")
        logging.info(f"The length of the data are: {len(spectro_intensity)}.")
        #logging.info(f"Intensity is: {spectro_intensity}.")
        logging.info(f"The type of data of y_dark is: {type(y_dark)}.")
        logging.info(f"the length of the data are: {len(y_dark)}.")
        #logging.info(f"Dark intensity is: {y_dark}.")
        logging.info(f"The type of data of intensity dar corrected is: {type(intensity_dark_corr)}.")
        logging.info(f'the length of the data are: {len(intensity_dark_corr)}.')
        #logging.info(f"Intensity dark corrected is: {intensity_dark_corr}.")

        self.send_data_to_user.emit("Correcting the power from dark.")
        """Correct the power from dark then calculate its square value."""
        power_dark_corr = (photodiode_power_mean) - (power_dark_mean)
        square_power_dark_corr = power_dark_corr * power_dark_corr

        logging.info("Power dark correction => ok !")

        self.send_data_to_user.emit("Calculate the area of the dark corrected emission.")
        """Calculate and store the area of the dark corrected emission"""

        logging.info(f"Sample informations: {self.sample_info.keys()}), sample: {sample}.")
        logging.info({self.sample_info[f"{sample}"]["Area correction factor"]})

        corr_dye_area = self.sample_info[f"{sample}"]["Area correction factor"]

        logging.info(f"Corr dye area= {corr_dye_area}.")

        self.send_data_to_user.emit(f"Corrected fluorescence is {round(corr_dye_area,1)}.")
        full_corrected_area = rdf()._Emission_Correction(corr_dye_area, intensity_dark_corr,
                                                         spectro_wavelength, self.ccd_correction_factor, emission_filter=self.current_emission_filter)
        logging.info("Emission_Correction => ok !")
        logging.info(f"P²= {square_power_dark_corr}.")
        logging.info(f'F= {full_corrected_area}.')

        try:
            fp=full_corrected_area/square_power_dark_corr
            self.send_data_to_user.emit(f"P² = {round(square_power_dark_corr,1)}. \nF = {round(full_corrected_area,1)}. \nF/P² = {round(fp,1)}.")

            logging.info(f'F/P² = {fp}.')

        except:

            logging.error(f"Can't calculate {full_corrected_area}/{square_power_dark_corr} !")

        self.signal_to_plot.emit(spectro_wavelength, spectro_intensity, square_power_dark_corr, full_corrected_area, i, true_power)
        if full_corrected_area <= 0:
            full_corrected_area = 1

            logging.warning("Full_corrected_area< 0. New value is set to 1.")

        """Save the measured data in the self.experimental dictionary."""
        if i == 0 and new_measure == False:
            self.experimental_data.update({
                sample:
                    {f"{wavelength}":
                         {i:
                              {'wavelength': spectro_wavelength,
                               'raw intensity': raw_fluorescence,
                               'raw power': photodiode_power,
                               'intensity - dark': intensity_dark_corr,
                               'raw power - dark': power_dark_corr,
                               'fully corrected area': full_corrected_area
                               },
                          "dark": y_dark
                          }
                     }
            })
        elif i != 0 and new_measure == False:
            self.experimental_data[f"{sample}"][f"{wavelength}"].update(
                {i:
                     {'wavelength': spectro_wavelength,
                      'raw intensity': raw_fluorescence,
                      'raw power': photodiode_power,
                      'intensity - dark': intensity_dark_corr,
                      'raw power - dark': power_dark_corr,
                      'fully corrected area': full_corrected_area
                      },
                 "dark": y_dark
                 }
            )
        elif new_measure == True:
            measure_number=self.data_for_new_measure[i]

            self.experimental_data[f"{sample}"][f"{wavelength}"].update(
                {measure_number:
                     {'wavelength': spectro_wavelength,
                      'raw intensity': raw_fluorescence,
                      'raw power': photodiode_power,
                      'intensity - dark': intensity_dark_corr,
                      'raw power - dark': power_dark_corr,
                      'fully corrected area': full_corrected_area
                      },
                 "dark": y_dark
                 }
            )

        logging.info(f"self.experimental for {sample} at {wavelength} and {true_power} mW updated.")

        self.send_data_to_user.emit(f"Experimental data of {sample} at {wavelength} nm and {round(float(true_power),5)} updated.")
        list_of_fully_corrected_area.append(full_corrected_area)
        list_of_square_dark_corr_power.append(square_power_dark_corr)

        logging.info(f"List of square dark corr power is: {list_of_square_dark_corr_power}.")
        logging.info(f"List of fully corrected area is: {list_of_fully_corrected_area}.")

        sleep(0.5)
        i += 1
        temporary_wavelength_list.append(spectro_wavelength)
        temporary_intensity_list.append(intensity_dark_corr)

        logging.info("Measure step 3 ended with success => ok !")

        return temporary_intensity_list, temporary_wavelength_list, i, list_of_fully_corrected_area, list_of_square_dark_corr_power, full_corrected_area

    def create_dictionary_processed_data(self, sample, wavelength, slope, coeff, log_square_fluo, log_square_power, fit,
                                         full_corrected_area, sigma2_phi, sigma2):
        """
        Save processed data in the dictionary self.processed_data
        :param sample: str
        :param wavelength: int
        :param slope: float
        :param coeff: float
        :param log_square_fluo: array
        :param log_square_power: array
        :param fit: float
        :param full_corrected_area: int
        :param sigma2_phi: float
        :param sigma2: float
        """
        self.processed_data[sample] = {
            f"{wavelength}":
                {
                    "quadraticity":
                        {
                            'slope': slope,
                            'coeff': coeff,
                            'log(F)': log_square_fluo,
                            'log(P²)': log_square_power,
                            'fit': fit
                        },
                    'full_corrected_area': full_corrected_area,
                    "S2F": sigma2_phi,
                    "S2": sigma2
                }
        }

    def update_dictionary_processed_data(self, sample, wavelength, slope, coeff, log_square_fluo, log_square_power, fit,
                                         full_corrected_area, sigma2_phi, sigma2):
        """
        Update the dictionary self.processed_data
        :param sample: str
        :param wavelength: int
        :param slope: float
        :param coeff: float
        :param log_square_fluo: array
        :param log_square_power: array
        :param fit: float
        :param full_corrected_area:int
        :param sigma2_phi: float
        :param sigma2: float

        """
        self.processed_data.update({
            sample:
                {
                    f"{wavelength}":
                        {
                            "quadraticity":
                                {
                                    'slope': slope,
                                    'coeff': coeff,
                                    'log(F)': log_square_fluo,
                                    'log(P²)': log_square_power,
                                    'fit': fit
                                },
                            'full_corrected_area': full_corrected_area,
                            "S2F": sigma2_phi,
                            "S2": sigma2
                        }
                }
        })

class MyApp(QWidget):
    """
    Main application for TPEF measurement.
    """

    def __init__(self):
        super().__init__()
        uic.loadUi('Enki.ui', self)
        self.setWindowTitle('Live Spectro')
        create_json()
        self.show_sample_info()
        create_samples("")
        self.refresh_data()
        self.processed_data = {}
        self.clean_labels()
        #self.label_Logo_NMP.setPixmap(QPixmap("Logo NMP.png"))
        self.folder_root= Path(get_path())
        save().create_angle_power_dict_for_measure()

        """ Page 1 - Initialization of the Calibration"""

        """Variables"""
        self.background = 0
        self.integration = self.integration_time_value
        self.intensity_int = 0
        self.power_int = 0
        self.ccd_correction_dye = None
        self.text = []
        self.user_info_calibration = []
        self.list_measuring_power = []
        self.list_measuring_time = []
        self.speed_of_angle_tuning = "fine"
        self.json_file_location = {}

        """Buttons"""
        self.btn_go.clicked.connect(self.run_loop_spectro)
        self.btn_stop.clicked.connect(self.stop_it)
        self.ckb_super_fast.clicked.connect(self.super_fast)
        self.btn_reset.clicked.connect(self.reset)
        self.btn_zero_power.clicked.connect(self.set_zero_power)
        self.btn_zero_intensity.clicked.connect(self.set_zero_intensity)
        self.btn_jod.clicked.connect(self.select_jod)
        self.btn_dmans.clicked.connect(self.select_dmans)
        self.btn_calibration.clicked.connect(self.run_ccd_correction)
        self.btn_set_wavelength.clicked.connect(self.wavelength_tuning)
        self.btn_experiment_status.setText("Safe")
        self.btn_experiment_status.setStyleSheet("background-color:green")
        self.btn_experiment_status_2.setText("Safe")
        self.btn_experiment_status_2.setStyleSheet("background-color:green")

        """List"""
        self.ls_terminal_calibration.clear()
        self.ls_terminal_calibration.setAutoScroll(True)

        """LCD Number"""
        self.lcd_power.display(self.power_int)
        self.lcd_intensity.display(self.intensity_int)

        """ProgressBar"""

        self.progressbar_calibration

        """ Page 2 - Initialization of the Parameters"""

        """Variables"""
        self.fluorescence_ref_dye = None
        self.ref_concentration = None
        self.checklist = {"calibration": False, "reference": False, "samples": True, "power": False,
                          "wavelength": False}
        self.user_info_parameters = []

        """Dictionnary"""
        self.sample_info = {}
        self.solvents = {"0": "Water", "1": "Methanol", "2": "Ethanol", "3": "DMSO", "4": "DMF", "5": "Acetonitrile",
                         "6": "Acetone", "7": "1,2-Dichloroethane", "8": "DCM", "9": "CHCl3", "10": "THF",
                         "11": "Acetic acid", "12": "1,4-Dioxane", "13": "Toluene", "14": "Cyclohexane",
                         "15": "Heptane"}
        self.refractive_index = {"Water": 1.33299, "Methanol": 1.32880, "Ethanol": 1.3611, "DMSO": 1.4770,
                                 "DMF": 1.4303, "Acetonitrile": 1.3441, "Acetone": 1.3588, "1,2-Dichloroethane": 1.4448,
                                 "DCM": 1.4242, "CHCl3": 1.44590, "THF": 1.4050, "Acetic acid": 1.3716,
                                 "1,4-Dioxane": 1.4224, "Toluene": 1.4961, "Cyclohexane": 1.4266, "Heptane": 1.3878}

        """Buttons"""
        self.btn_set_fluorescein.clicked.connect(self.select_fluorescein)
        self.btn_set_nile_red.clicked.connect(self.select_nile_red)
        self.btn_create_sample_number.clicked.connect(self.create_sample)
        self.btn_delete_sample.clicked.connect(self.remove_sample)
        self.btn_sample_info.clicked.connect(self.show_sample_info)
        self.btn_open.clicked.connect(self.open_emission_root)

        """Check Buttons"""
        self.btn_check_calibration.setChecked(False)
        self.btn_check_reference.setChecked(False)
        self.btn_check_samples.setChecked(False)
        self.btn_check_wavelength.setChecked(False)
        self.btn_check_power.setChecked(False)
        self.btn_check_opo.clicked.connect(self.set_opo)

        """List"""
        self.ls_list_of_samples.clear()
        self.ls_terminal_parameters.clear()

        """Misc"""
        self.show_solvent()

        """spinbox"""
        self.spb_set_wavelength_min.valueChanged.connect(self.set_wavelength)
        self.spb_set_wavelength_max.valueChanged.connect(self.set_wavelength)
        self.spb_set_wavelength_pitch.valueChanged.connect(self.set_wavelength)
        self.spb_set_power_min.valueChanged.connect(self.set_power)
        self.spb_set_power_max.valueChanged.connect(self.set_power)
        self.spb_set_number_of_power_measure.valueChanged.connect(self.set_power)

        """QDial"""
        self.dial_power.valueChanged.connect(lambda:self.set_power_laser())

        """Line edit"""
        self.le_ref_concentration.textChanged.connect(self.concentration)

        """ Page 3 - Initialization of the Measurement"""

        """Variables"""
        self.user_info_measurement = []
        self.order_of_measure = []
        self.sigma2_dict = {}
        self.wavelength_list = []
        self.sigma2_list = []
        self.sigma2phi_list = []
        self.experimental_data = {}
        self.processed_data = {}
        self.counter_results = {}
        self.root_dict = {}
        self.done_wavelengths = []
        self.first_measure= True


        """Intruments"""
        self.spectro= talk_to_spectro

        """Timer"""
        self.timer_spectro = QTimer()
        self.timer_powermeter = QTimer()
        # self.timer.setInterval(30) <= minimal value to allow stop function working
        self.timer_spectro.setInterval(100)
        self.timer_powermeter.setInterval(30)
        self.timer_spectro.timeout.connect(self.update_graphs_spectro)
        self.timer_powermeter.timeout.connect(self.update_graphs_power)
        self.faster = False

        """Button"""
        self.btn_new_measure_data.clicked.connect(self.new_measure_of_data)
        self.btn_remove_data.clicked.connect(self.remove_data)
        self.btn_refresh.clicked.connect(self.refresh_data)
        self.btn_set_sample_wavelength.clicked.connect(self.show_data_of_sample)
        self.btn_run_tpa.clicked.connect(self.run_tpa_measurements)
        self.btn_show_sigma_2.clicked.connect(self.plot_s2_from_signal)
        self.btn_simulation.clicked.connect(self.simulation)

        """List"""
        self.ls_terminal_information_for_users_measurement.clear()
        self.ls_terminal_information_for_users_measurement.setAutoScroll(True)
        """ Page 4 - Advanced Calibration"""

        """Button"""

        self.btn_filter_FF01650.clicked.connect(self.select_filter_FF01650)
        self.btn_filter_FF01650.setChecked(True)
        self.select_filter_FF01650()
        self.btn_filter_E650sp2p.clicked.connect(self.select_filter_E650sp2p)
        self.btn_filter_E750sp2p.clicked.connect(self.select_filter_E750sp2p)
        self.btn_filter_E650.clicked.connect(self.select_filter_E650)
        self.btn_filter_FES700.clicked.connect(self.select_filter_FES700)
        self.btn_filter_FES750.clicked.connect(self.select_filter_FES750)
        self.btn_filter_FES800.clicked.connect(self.select_filter_FES800)

        """Variables"""

        self.current_emission_filter = "FF01-650SP25"

        """Properties"""
    @property
    def power_max(self):
        return self.spb_set_power_max.value()
    @property
    def power_min(self):
        return self.spb_set_power_min.value()
    @property
    def power_pitch(self):
        return self.spb_set_number_of_power_measure.value()
    @property
    def scan_value(self):
        return self.spb_scan.value()
    @property
    def integration_time_value(self):
        return self.spb_integration.value()
    @property
    def wavelength_tuning_value(self):
        return self.spb_wavelength_tuning.value()
    @property
    def wavelength_min(self):
        return self.spb_set_wavelength_min.value()
    @property
    def wavelength_max(self):
        return self.spb_set_wavelength_max.value()
    @property
    def wavelength_pitch(self):
        return self.spb_set_wavelength_pitch.value()
    """
    Methods for all
    """
    def communication_with_user(self, data_to_show, clear=False, calibration=False, information=False, measurement=False):
        """
        Send data to user using the QListwidgets ls_terminal_calibration, ls_terminal_information_for_users_parameters
        or ls_terminal_information_for_users_measurement
        :param data_to_show: str
        :param clear: bool
        :param calibration: bool
        :param information: bool
        :param measurement: bool
        """
        if calibration :
            if clear :
                self.ls_terminal_calibration.clear()
            self.ls_terminal_calibration.addItem(str(data_to_show))
            self.ls_terminal_calibration.update()
        elif information :
            if clear:
                self.ls_terminal_information_for_users_parameters.clear()
            self.ls_terminal_information_for_users_parameters.addItem(str(data_to_show))
            self.ls_terminal_information_for_users_parameters.update()
        elif measurement :
            if clear:
                self.ls_terminal_information_for_users_measurement.clear()
            self.ls_terminal_information_for_users_measurement.addItem(str(data_to_show))
            self.ls_terminal_information_for_users_measurement.update()

        logging.info("Send data to user interface => ok !")

    """
    Calibration
    """

    def set_power_laser(self):
        value = self.dial_power.value()
        self.label_dial_power.setText(f"Required angle value= {value}°.")
        print(value)
        tuner.RotationMount().spin_to_position(int(value))
        time.sleep(0.2)
        angle=tuner.RotationMount().get_angle()
        print(angle)
        self.label_dial_power.setText(f"Real angle value= {angle}°.")

    def all(self):
        """
        Send data and update data to PlotWidgets Spectre & Power
        """
        self.update_graphs_spectro()
        self.update_graphs_power()

        logging.info("Spectre & Power plotwidgets updated => ok !")

    def close_photodiode(self):
        """
        Close the photodiode
        """
        if SIMULATION_FOR_DEBUG == False:
            with nidaqmx.Task() as task:
                task.ai_channels.add_ai_voltage_chan("Dev1/ai1")
                task.close()
        else:
            print("Close")

        logging.info("Photodiode closed.")

    def get_ccd_corr_dye(self):
        """
        Collect and return data about the dye used for CCD correction
        :return: int (excitation wavelength), int (power of the laser), str (color for plot)
        """
        ccd_ref_dict = References_data.ccd_corr_dye()
        try :
            for ref_name in ccd_ref_dict:
                if self.ccd_correction_dye == str(ref_name):
                    wavelength = ccd_ref_dict[f"{ref_name}"]["wavelength"]
                    power_for_ccd_dye = ccd_ref_dict[f"{ref_name}"]["power"]
                    color = ccd_ref_dict[f"{ref_name}"]["color"]

                    logging.info(f"The selected reference is {ref_name}.")
                    logging.info(f"Wavelength= {wavelength} nm, type: {type(wavelength)}.")
                    logging.info(f"power= {power_for_ccd_dye} mW, type: {type(power_for_ccd_dye)}")

            return wavelength, power_for_ccd_dye, color
        except UnboundLocalError:

            logging.error("Can't get informations about selected CCD correction dye !")

            return None, None, "g"
    def graphs(self, x, y):
        """
        Plot emission spectrum on PlotWidget Spectre.
        :param x: array (wavelength)
        :param y: array (intensity)
        """
        logging.info("Plotinf data on Spectre.")

        x = np.array(x)
        y = np.array(y)
        self.Spectre.clear()
        self.Spectre.addLegend()
        wavelength, power, color=self.get_ccd_corr_dye()
        self.Spectre.setLabel('left', units="Fluorescence intensity / cps")
        self.Spectre.setLabel('bottom', units="Wavelength / nm")
        #self.Spectre.setYRange(0, 16000)
        self.Spectre.setXRange(320, 1100)
        try :
            self.Spectre.plot(x, y, pen= f"{color}",name= f"Intensity of {self.ccd_correction_dye}.", show=True)

        except UnboundLocalError:
            self.Spectre.plot(x, y, pen=f"{color}", name="Intensity", show=True)

            logging.warning("Invalid self.ccd_correction_dye for plotting !")

        self.intensity_integration(x, y)
    def graphs_p(self, time, power):
        """
        Plot laser power stability on PlotWidget Power.
        :param time: array
        :param power: array
        """

        logging.info("Plotting data on Power.")

        self.list_measuring_power.append(power)
        self.list_measuring_time.append(time)
        self.Power.setLabel("left", units= "power / mV")
        self.Power.setLabel("bottom", units="Time / s")
        self.Power.plot(self.list_measuring_time, self.list_measuring_power, pen= "r", symbol= "o", symbolPen= "r", symbolBrush= 0.5, show=True)
        self.power_integration(power)
    def intensity_integration(self, x, y):
        """
        Send to the user using the lcd_intensity label the value of the fluorescence integral.
        :param x: array (wavelength)
        :param y: array (intensity)
        """
        self.intensity_int = np.trapz(y, x)
        self.intensity_int /= 1e5
        self.lcd_intensity.display(self.intensity_int - self.background)

        logging.info("Calculation of the integral of intensity => ok !")

    def power_integration(self, y):
        """
        Send to the user using the lcd_power label the value of the power
        :param y: float
        :return:
        """
        self.power_int = y #np.mean(y)
        self.lcd_power.display(self.power_int * 1000)
    def read_photodiode(self):
        """
        Read and return data from the photodiode
        :return: float
        """
        if SIMULATION_FOR_DEBUG == False:
            with nidaqmx.Task() as task:
                task.ai_channels.add_ai_voltage_chan("Dev1/ai1", max_val=5, min_val=-5)
                #data = task.read()
                in_stream = task.in_stream

                data = in_stream.read(number_of_samples_per_channel=self.integration_time_value)
                print('1 Channel N Samples Read: ')
                #data = task.read(number_of_samples_per_channel=10)
                print(len(data), type(data))
                print(data)
                data=np.mean(data)
                print(data)


        else:
            data = np.random.randint(10)
        return data
    def reset(self):
        """
        Clear ls_terminal_calibration, PlotWiget Power and lists of measured time & power

        """
        self.list_measuring_power.clear()
        self.list_measuring_time.clear()
        self.Power.clear()
        self.ls_terminal_calibration.clear()

        logging.info("Power calibration interface reset by user.")

    def run(self, intensity_dark):
        """
        Read the spectro data and correct the intensity from the dark noise.
        :param intensity_dark: array (intensity dark noise)
        :return: array (wavelength), array (intensity with correction of the dark noise)
        """
        logging.info("Running spectro acquisition.")

        x, y , raw_fluo= self.spectro.Maya(int(self.integration_time_value))._data_read_scan(number_of_scan=1, intensity_dark=intensity_dark)

        logging.info("Collecting data from spectro => ok !")
        logging.info(f"Wavelength= {x}, intensity= {y}.")

        return x, y
    def run_ccd_correction(self):
        """
        Perform the measurement and the calculation of the correction of the USB spectro using a fluorescent
        dye. Set the button correction check to True
        :return: True
        """
        logging.info("Starting CCD correction.")

        self.btn_experiment_status.setText("Measure in progress...")
        self.btn_experiment_status.setStyleSheet("background-color:red")
        self.progressbar_calibration.setValue(0)
        self.ls_terminal_calibration.clear()

        """Check that the user has selected a CCD reference for the correction."""
        if self.ccd_correction_dye == None:
            self.ls_terminal_calibration.clear()
            self.communication_with_user(data_to_show="Please select a dye for the CCD correction", calibration=True)

            logging.error("No CCD dye selected.")

            return False
        else:
            """Initialization of the spectro, setting the integration time, check that the shutter of the laser
            is closed, otherwise close it."""

            logging.info(f"Fluorescence of {self.ccd_correction_dye} is measured.")

            self.communication_with_user(data_to_show=f"Starting measurement of {self.ccd_correction_dye}.", calibration=True)
            shutter_state = fs_laser().queryShutterStatus()
            if shutter_state == b's=1':
                fs_laser().closeShutterBlocking()

            logging.info("Laser shutter is close.")

            self.progressbar_calibration.setValue(10)
            sleep(0.5)
            wavelength, power_for_ccd_dye, color = self.get_ccd_corr_dye()
            self.progressbar_calibration.setValue(20)

            """Send the suitable wavelength to the laser, determine the angle of the half-wave plate to get a power
            ~ 100 mW, then set the half-wave plate to this angle."""
            fs_laser().setWavelengthBlocking(wavelength)

            logging.info(f"Laser tuned to {wavelength} nm")

            self.communication_with_user(data_to_show=(f"Tuning the laser to {wavelength} nm."), calibration=True)
            self.progressbar_calibration.setValue(30)
            sleep(0.5)
            angle_min, angle_max, all_data = tuner.RotationMount().determine_angle_for_power(power_for_ccd_dye, 90, 0, wavelength)

            logging.info(f"Angle_max is: {angle_max}.")

            tuner.RotationMount().spin_to_position(angle_max)

            self.communication_with_user(data_to_show=f"The angle of halfwave is {angle_max}.", calibration=True)
            sleep(0.5)
            self.progressbar_calibration.setValue(40)

            """Dark spectrum for the CCD correction"""

            logging.info("Measuring the dark.")

            self.communication_with_user(data_to_show="Measuring the dark..", calibration=True)
            sleep(0.5)
            x_dark, y_dark, raw_dark_intensity = self.spectro.Maya(int(self.integration_time_value))._data_read_scan(number_of_scan=1, intensity_dark=0)

            logging.info("Collecting data from spectro => ok !")

            self.graphs(x_dark, y_dark)
            self.progressbar_calibration.setValue(50)
            self.communication_with_user(data_to_show="Dark measurement successful.", calibration=True)
            """Fluorescence intensity of the ref for the CCD correction"""
            self.communication_with_user(data_to_show=("Waiting for laser stabilization.. \n Opening the laser shutter."), calibration=True)
            fs_laser().openShutterBlocking()

            logging.info("Laser shutter is open.")

            self.progressbar_calibration.setValue(60)
            time_shift = 0
            for unit_time in range(15):
                sleep(1)
                self.progressbar_calibration.setValue(60 + time_shift)
                time_shift += 1
            self.communication_with_user(data_to_show=(f"Measuring the {self.ccd_correction_dye}.."), calibration=True)
            sleep(0.5)

            logging.info(f"Measuring fluorescence of {self.ccd_correction_dye}.")

            x, y_corr,  raw_dark_intensity = self.spectro.Maya(int(self.integration_time_value))._data_read_scan(number_of_scan=self.scan_value, intensity_dark=y_dark)

            logging.info("Collecting data from spectro => ok !")

            self.communication_with_user(data_to_show=(f"Measurement of {self.ccd_correction_dye} successful. \n Closing the laser shutter."), calibration=True)
            fs_laser().closeShutterBlocking()

            logging.info("Laser shutter is closed.")

            self.graphs(x, y_corr)
            self.progressbar_calibration.setValue(80)

            """Calculation of the CCD correction factor and send a message to the user and display the spectrum."""
            #y_corr = y - y_dark

            self.graphs(x, y_corr)
            wavelength, self.ccd_correction_factor = rdf()._Spectro_Correction_Factor(x, y_corr, emission_filter=self.current_emission_filter, dye_for_correction=self.ccd_correction_dye)
            self.communication_with_user(data_to_show=("Calculation of the correction factor."), calibration=True)
            self.progressbar_calibration.setValue(85)
            if type(self.ccd_correction_factor) == bool:
                self.communication_with_user(data_to_show=f"Error with the calibration step ! \n Please check that \n the emission filter select \n in advanced calibration fits \n with {self.ccd_correction_dye} !", calibration=True)

                logging.error("Error with the calibration step ! Check the emission filter selected !")
                self.progressbar_calibration.setValue(100)
                self.btn_experiment_status.setText("Safe")
                self.btn_experiment_status.setStyleSheet("background-color:green")
                return False
            self.btn_check_calibration.setChecked(True)
            self.checklist["calibration"] = True
            self.progressbar_calibration.setValue(90)
            self.communication_with_user(data_to_show=f"Calibration with {self.ccd_correction_dye} successful !", calibration=True)

            #logging.info(f"CCD _correction factor is: {self.ccd_correction_factor}.")

            save().figure(x=x, y=y_corr, name=self.ccd_correction_dye, sample="CCD correction dye")
            save().figure(x=wavelength, y=self.ccd_correction_factor, name= "CCD correction curve", sample="CCD correction dye")
            self.progressbar_calibration.setValue(100)
            self.btn_experiment_status.setText("Safe")
            self.btn_experiment_status.setStyleSheet("background-color:green")
            #self.close_photodiode() #remove if issues
            return True
    def run_loop_spectro(self):
        """
        First measure the dark noise of the spectro then measure the power and read the spectro for real time.
        """
        self.btn_experiment_status.setText("Measure in progress...")
        self.btn_experiment_status.setStyleSheet("background-color:red")
        self.Power.clear()
        talk_to_spectro.Maya(int(self.integration_time_value))

        logging.info(f"Integration time= {self.integration_time_value} ms.")

        logging.info("Measuring dark for live acquisition.")

        self.x_dark_live, self.y_dark_live = self.run(intensity_dark=0)
        sleep(0.5)
        if self.faster == False:
            self.timer_spectro.setInterval(int(self.integration_time_value))
        else:
            self.timer_spectro.setInterval(30)
        self.ls_terminal_calibration.addItem("Connection to USB Spectro successful")
        self.ls_terminal_calibration.addItem(f"Integration time set to {self.integration_time_value} ms")
        self.ls_terminal_calibration.addItem(30 * "#")
        fs_laser().openShutterBlocking()
        sleep(0.1)
        self.start_measure_power = chrono().run()
        self.timer_powermeter.start()
        self.timer_spectro.start()
    def select_dmans(self):
        """
        Makes the CCD correction dye = DMANs
        """
        self.ccd_correction_dye = "DMANs"
        self.ls_terminal_calibration.addItem("The dye selected is DMANs")

        logging.info("DMANS was selected as CCD dye.")

    def select_jod(self):
        """
        Makes the CCD correction dye = JOD
        """
        self.ccd_correction_dye = "JOD"
        self.ls_terminal_calibration.addItem("The dye selected is JOD")

        logging.info("JOD was selected as CCD dye.")

    def set_zero_intensity(self):
        """
        Set the background value of the power to 0.
        """
        self.background = self.intensity_int
        self.lcd_intensity.display(self.intensity_int - self.background)
        self.ls_terminal_calibration.addItem(30 * "#")
        self.ls_terminal_calibration.addItem("The zero is done !")
    def set_zero_power(self):
        """
        Send 0 to module power_integration
        """
        self.power_integration(0)
        self.ls_terminal_calibration.addItem(30 * "#")
        self.ls_terminal_calibration.addItem("The zero is done !")
    def stop_it(self):
        """
        Stop the acquisition of the spectro & photodiode data and close the shutter of the laser.
        """
        fs_laser().closeShutterBlocking()

        logging.info("Laser shutter is closed")

        self.timer_powermeter.stop()
        self.timer_spectro.stop()
        self.ls_terminal_calibration.addItem(30 * "#")
        self.ls_terminal_calibration.addItem("Pause")
        self.ls_terminal_calibration.addItem(30 * "#")
        self.list_measuring_power.clear()
        self.list_measuring_time.clear()
        self.btn_experiment_status.setText("Safe")
        self.btn_experiment_status.setStyleSheet("background-color:green")
        #self.close_photodiode()
    def super_fast(self):
        """
        Turns self.faster to True to set integration time to 30 ms.
        """
        self.faster = True

        logging.info("User selected super fast acquisition (30 ms).")

        return self.faster
    def update_graphs_power(self):
        """
        Read the photodiode and update the PlotWidget Power with the new value.
        """
        y = self.read_photodiode()
        x = chrono().get_time(self.start_measure_power)

        logging.info(f"updating Power, time= {x}, power= {y}.")

        self.graphs_p(x, y)
    def update_graphs_spectro(self):
        """
        Read spectro data and update the PlotWiget Spectre with the new datas
        """
        x, y = self.run(intensity_dark=self.y_dark_live)
        y_corr = y#-self.y_dark_live
        self.graphs(x, y_corr)

        logging.info(f"updating Spectre, wavelength= {x}, Fluorescence= {y_corr}.")

    def wavelength_tuning(self):
        """
        Tune the laser to the desire wavelength
        """
        wavelength = self.wavelength_tuning_value
        if wavelength > 679 and wavelength < 1081:
            fs_laser().setWavelengthBlocking(int(wavelength))

            logging.info(f"Laser tuned to {wavelength} nm.")

        else:

            logging.warning(f"The selected wavelength ({wavelength} nm) is not suitable for the laser.")

            return

    """
    Parameters
    """


    def add_sample(self):
        """
        Add a new sample to widget ls_list_of_samples
        """
        sample_name = self.le_set_sample_name.text()
        if sample_name == "":

            logging.warning("User do not enter name for the sample name.")

            return False
        else:
            sample = sample_data(sample_name=sample_name)
            result = sample.add_to_sample_name()
            if result:
                lw_item = QListWidgetItem(sample.sample_name)
                lw_item.setData(0x0100, sample)  # 0x0100 => QtCore.Qt.UserRole
                self.ls_list_of_samples.addItem(lw_item)
                self.le_set_sample_name.setText("")

        logging.info("A new sample was added.")

    def concentration(self):
        """
        Check that the user send a valid concentration.
        """
        try:
            self.ref_concentration = self.le_ref_concentration.text()

            if self.ref_concentration == "" or self.ref_concentration.replace(".", "", 1).replace(",", "", 1).replace(
                    "e", "", 1).isdigit() == False:
                self.user_info_parameters.clear()
                self.user_info_parameters.append("Concentration input is wrong")

                logging.warning("Concentration input is wrong !")

                self.ls_terminal_parameters.clear()
                for _ in self.user_info_parameters:
                    self.ls_terminal_parameters.addItem(_)
                return False

            if self.fluorescence_ref_dye != None:
                self.btn_check_reference.setChecked(True)
                self.checklist["reference"] = True
                self.ls_terminal_parameters.addItem(
                    f"The concentration of {self.fluorescence_ref_dye} is {self.ref_concentration} mol/L")

                logging.info(f"The concentration of {self.fluorescence_ref_dye} is {self.ref_concentration} mol/L")

        except ValueError:

            logging.error("Concentration_Value ERROR.")
            return

    def create_sample(self):
        """
        Create a new sample using the data send by the user : name, concentration, solvent, root of emission spectrum
        & phi. Save these data to the dictionary sample_info
        """
        solvent = None
        sample_name = self.le_set_sample_name.text()
        if sample_name == "None":
            self.user_info_parameters.clear()
            self.user_info_parameters.append("None cannot be used as sample name. \n Please change sample name.")

            logging.warning("None cannot be used as sample name, please change sample name !")

            self.ls_terminal_parameters.clear()
            for _ in self.user_info_parameters:
                self.ls_terminal_parameters.addItem(_)
            return False
        concentration = self.le_set_sample_concentration.text()
        if concentration == "" or concentration.replace(".", "", 1).replace(",", "", 1).replace("e", "",
                                                                                                1).isdigit() == False:
            self.user_info_parameters.clear()
            self.user_info_parameters.append("Concentration input is wrong")

            logging.warning("Concentration input is wrong !")

            self.ls_terminal_parameters.clear()
            for _ in self.user_info_parameters:
                self.ls_terminal_parameters.addItem(_)
            return False
        emission_root = self.sample_root[0]
        phi = self.le_set_fluo_qy.text()
        if phi == "" or phi.replace(".", "", 1).replace(",", "", 1).isdigit() == False:
            self.user_info_parameters.clear()
            self.user_info_parameters.append("Fluorescence quantum yield input is wrong")

            logging.warning("Fluorescence quantum yield input is wrong !")

            self.ls_terminal_parameters.clear()
            for _ in self.user_info_parameters:
                self.ls_terminal_parameters.addItem(_)
            return False
        elif float(phi) > 0 and float(phi) < 1:
            pass
        else:
            self.user_info_parameters.clear()
            self.user_info_parameters.append("Fluorescence quantum yield input is wrong")

            logging.warning("Fluorescence quantum yield input is wrong !")

            self.ls_terminal_parameters.clear()
            for _ in self.user_info_parameters:
                self.ls_terminal_parameters.addItem(_)
            return False
        for selected_solvent in self.ls_solvent.selectedItems():
            solvent = selected_solvent.data(0x0100)
        if solvent == None:
            self.user_info_parameters.clear()
            self.user_info_parameters.append("Please select a solvent")

            logging.warning("No solvent selected !")

            self.ls_terminal_parameters.clear()
            for _ in self.user_info_parameters:
                self.ls_terminal_parameters.addItem(_)
            return False
        emission_wavelength = []
        emission_intensity = []
        with open(emission_root, "r") as item:
            values = item.readlines()
            for item in values:
                t = item.replace("\t\n", "").split("\t")
                emission_wavelength.append(float(t[0]))  # first index of the list is wavelength
                emission_intensity.append(float(t[1]))  # second index of the list is fluorescence intensity
        corr_dye_area = rdf()._Read_Dyes_Emission_Spectrum(root=emission_root, sample_name=sample_name, emission_filter=self.current_emission_filter)
        self.counter_results[sample_name] = True

        logging.info(f"Corr_dye_area: {corr_dye_area}.")

        self.add_sample()
        self.sample_info[f"{sample_name}"] = {"emission": f"{emission_root}", "concentration": f"{concentration}",
                                              "solvent": f"{solvent}", "phi": f"{phi}",
                                              "emission_wavelength": emission_wavelength,
                                              "emission_intensity": emission_intensity,
                                              "Area correction factor": corr_dye_area,
                                              "CCD reference": self.ccd_correction_dye, "state": True}
        self.le_set_sample_concentration.setText("")
        self.le_set_fluo_qy.setText("")
        self.btn_check_samples.setChecked(True)
        list_of_sample = [self.fluorescence_ref_dye]
        for key in self.sample_info:
            list_of_sample.append(key)
        list_of_sample.append(self.fluorescence_ref_dye)
        create_samples(list_of_samples= list_of_sample)
        self.refresh_data()
    def open_emission_root(self):
        """
        Get the path of the emission spectrum file and store it in self.sample_root
        """
        self.sample_root = QFileDialog.getOpenFileName()

        logging.info(f"{self.sample_root[0]}.")

    def remove_sample(self):
        """
        Remove the selected sample in the ls_list_of_samples.
        """
        for selected_item in self.ls_list_of_samples.selectedItems():
            sample_to_remove = selected_item.data(0x0100)
            sample_data(sample_name=f"{sample_to_remove}").remove_sample()
            del self.sample_info[f"{sample_to_remove}"]
            sample_to_remove.remove_sample()
            self.ls_list_of_samples.takeItem(self.ls_list_of_samples.row(selected_item))
        self.refresh_data()

        logging.info(f"The sample {selected_item} was removed from sample list..")

    def select_fluorescein(self):
        """
        Select fluorescein as fluorescent dye for TPEF calculation
        :return: True
        """
        self.fluorescence_ref_dye = "fluorescein"
        self.ls_terminal_parameters.addItem("The ref dye is fluorescein")

        logging.info("Fluorescein was selected.")

        return True
    def select_nile_red(self):
        """
        Select Nile red as fluorescent dye for TPEF calculation
        :return: True
        """
        self.fluorescence_ref_dye = "NR"
        self.ls_terminal_parameters.addItem("The ref dye is nile red")

        logging.info("Nile Red was selected.")

        return True
    def set_power(self):
        """
        Check the checkbox power
        """
        self.btn_check_power.setChecked(True)
        self.checklist["power"] = True

        logging.info("btn_check_power was checked.")

    def set_opo(self):
        """
        Allow to extend wavelength range to 1400 is user check the OPO checkbox
        """
        if self.btn_check_opo.isChecked():

            logging.info("User selected OPO.")

            self.spb_set_wavelength_min.setMaximum(1400)
            self.spb_set_wavelength_max.setMaximum(1400)

        else:
            print("OPO uncheck")
            self.spb_set_wavelength_min.setMaximum(1080)
            self.spb_set_wavelength_max.setMaximum(1080)
    def set_wavelength(self):
        """
        Check the checkbox wavelength.
        """
        self.wavelength_range = cfg().wavelength(self.wavelength_min, self.wavelength_pitch, self.wavelength_max)

        self.refresh_wavelength(self.wavelength_range)
        self.btn_check_wavelength.setChecked(True)

        self.checklist["wavelength"] = True
    def show_sample_info(self):
        """
        Send to the user using the ls_terminal_parameters, data about selected sample
        """
        for selected_samples in self.ls_list_of_samples.selectedItems():
            sample = selected_samples.data(0x0100)
            concentration = self.sample_info[f"{sample}"]["concentration"]
            phi = self.sample_info[f"{sample}"]["phi"]
            solvent = self.sample_info[f"{sample}"]["solvent"]
            self.user_info_parameters.clear()
            self.user_info_parameters.append(
                f"Sample name : {sample}\n Concentration = {concentration}\n Fl. Quantum yield = {phi}\n Solvent = {solvent}\n")
            self.ls_terminal_parameters.clear()
            for _ in self.user_info_parameters:
                self.ls_terminal_parameters.addItem(_)
            self.user_info_parameters.clear()

    def show_solvent(self):
        """
        Set the entries in the spinbox solvent from dictionary self.solvents.
        """
        for key, value in self.solvents.items():
            lw_item = QListWidgetItem(self.solvents[key])
            lw_item.setData(0x0100, self.solvents[key])  # 0x0100 => QtCore.Qt.UserRole
            self.ls_solvent.addItem(lw_item)

    """
    Measurement _ Measure
    """
    def create_process_data_files(self):
        """
        Create the process data file for each sample.
        """
        for key in self.sample_info:
            self.experimental_data[key] = {}
            self.processed_data[key] = {}
            save().write_only(dictionary_process= self.processed_data[key], sample= key)

            logging.info("Process data file was successfully created.")

    def disable_parameters_widgets(self):
        """
        Disable all key widgets to avoid user to manipulate data during measurements.
        :return: True
        """

        logging.info("Parameters widget are disabled.")

        """Buttons"""
        """Calibration"""
        self.btn_jod.setEnabled(False)
        self.btn_dmans.setEnabled(False)
        self.btn_zero_power.setEnabled(False)
        self.btn_zero_intensity.setEnabled(False)
        self.btn_stop.setEnabled(False)
        self.btn_go.setEnabled(False)
        self.btn_reset.setEnabled(False)
        self.btn_calibration.setEnabled(False)
        self.btn_set_wavelength.setEnabled(False)

        """Parameters"""
        self.btn_set_fluorescein.setEnabled(False)
        self.btn_set_nile_red.setEnabled(False)
        self.btn_create_sample_number.setEnabled(False)
        self.btn_delete_sample.setEnabled(False)
        self.btn_sample_info.setEnabled(False)
        self.btn_open.setEnabled(False)
        self.btn_check_calibration.setEnabled(False)
        self.btn_check_reference.setEnabled(False)
        self.btn_check_samples.setEnabled(False)
        self.btn_check_power.setEnabled(False)
        self.btn_check_wavelength.setEnabled(False)
        self.btn_new_measure_data.setEnabled(False)
        self.btn_remove_data.setEnabled(False)
        self.btn_run_tpa.setEnabled(False)

        """Measurement"""
        self.btn_reset_power_angle_conversion_file.setEnabled(False)

        """Check box"""
        """Calibration"""
        self.ckb_super_fast.setEnabled(False)

        """Spin box"""
        """Calibration"""
        self.spb_integration.setEnabled(False)
        self.spb_wavelength_tuning.setEnabled(False)

        """Parameters"""
        self.spb_set_wavelength_min.setEnabled(False)
        self.spb_set_wavelength_max.setEnabled(False)
        self.spb_set_wavelength_pitch.setEnabled(False)
        self.spb_set_power_min.setEnabled(False)
        self.spb_set_power_max.setEnabled(False)
        self.spb_set_number_of_power_measure.setEnabled(False)

        """Line edit"""
        """Parameters"""
        self.le_ref_concentration.setEnabled(False)
        self.le_set_sample_name.setEnabled(False)
        self.le_set_sample_concentration.setEnabled(False)
        self.le_set_fluo_qy.setEnabled(False)

        """List"""
        self.ls_solvent.setEnabled(False)

        return True
    def enable_widget(self):
        """
        Liberate the widgets.
        """

        logging.info("Parameters widget are enabled.")

        """Buttons"""
        """Calibration"""

        self.btn_go.setEnabled(True)
        self.btn_stop.setEnabled(True)
        self.btn_reset.setEnabled(True)
        self.btn_set_wavelength.setEnabled(True)
        self.btn_zero_intensity.setEnabled(True)
        self.btn_zero_power.setEnabled(True)

        """Measurement"""
        self.btn_new_measure_data.setEnabled(True)
        self.btn_remove_data.setEnabled(True)

        """Spin box"""
        """Calibration"""
        self.spb_wavelength_tuning.setEnabled(True)
        self.spb_integration.setEnabled(True)

        """Check box"""
        """Calibration"""
        self.ckb_super_fast.setEnabled(True)
    def enable_widgets_after_run(self):
        """
        Liberate widgets after run.
        """

        logging.info("After run widgets enabled.")

        self.btn_set_wavelength.setEnabled(True)
        self.ckb_super_fast.setEnabled(True)
        self.spb_integration.setEnabled(True)
        self.spb_wavelength_tuning.setEnabled(True)
        self.btn_new_measure_data.setEnabled(True)
        self.btn_remove_data.setEnabled(True)
    def measurement(self, sample, wavelength, new_measure, data_for_new_measure,opo=False):
        """
        Manage the start and ends of the measurement.
        :param sample: str, sample name
        :param wavelength: int, excitation wavelength
        :param new_measure: bool, measuring again a wavelength/sample
        :param data_for_new_measure: dict, old data if new_measure == True
        """
        self.ls_measure_data.clear()

        logging.info("Starting measurement of TPA.")

        if sample == "None":
            print("Please select a sample")
            self.enable_widgets_after_run()
            self.btn_run_tpa.setEnabled(True)

            logging.warning("No sample selected !")

            return
        else:
            self.btn_experiment_status_2.setText("Measure in progress...")
            self.btn_experiment_status_2.setStyleSheet("background-color:red")
            self.thread = QThread()
            self.worker = Measure(sample=sample,
                                  wavelength=wavelength,
                                  integration_time=self.integration_time_value,
                                  dictionary_of_sample_info=self.sample_info,
                                  dictionary_of_sample_experimental_data=self.experimental_data,
                                  dictionary_of_sample_processed_data=self.processed_data,
                                  fluorescence_ref_dye=self.fluorescence_ref_dye,
                                  dictionary_of_sample_root=self.root_dict,
                                  ccd_correction_dye=self.ccd_correction_dye,
                                  ccd_correction_factor=self.ccd_correction_factor,
                                  maximum_power_requested=self.power_max,
                                  minimum_power_requested=self.power_min,
                                  power_pitch=self.power_pitch,
                                  number_of_scan=self.scan_value,
                                  new_measure=new_measure,
                                  data_for_new_measure=data_for_new_measure,
                                  current_emission_filter=self.current_emission_filter,
                                  opo_measure=self.opo_measure_state,
                                  real_opo_wavelength=self.opo_wavelength
                                  )

            self.worker.moveToThread(self.thread)
            self.worker.signal_to_plot.connect(self.plot_from_signal)
            self.worker.enable_widgets.connect(self.enable_widgets_after_run)
            self.worker.signal_to_pop_up.connect(self.pop_up_end_of_measure)
            self.worker.save_quad.connect(self.save_quad_from_signal)
            self.worker.save_fluo_vs_power.connect(self.save_fluo_vs_power)
            self.worker.clear_plot.connect(self.clear_plot_from_signal)
            self.worker.send_data_to_user.connect(self.show_data_to_user)
            self.worker.plot_to_save.connect(self.save_plot_from_signal)
            self.worker.progress_bar_status.connect(self.change_progress_bar_measurment_status)
            self.worker.measure_finished.connect(self.measure_end)
            print("Thread ready !")
            self.thread.started.connect(self.worker.run_measure)
            self.worker.finished.connect(self.thread.quit)
            print("Thread ready step 2 !")
            self.thread.start()
            print("Thread ready step 3 !")
    def change_progress_bar_measurment_status(self, value):
        """
        Set the value of the progressbar_measurement
        :param value: int, range 0-100
        :return: True
        """
        self.progressbar_measurement.setValue(int(value))
        return True
    def measure_end(self, sample_info, experimental_data, processed_data, root_dict):
        """
        Store updated dictionary in the gui after measurement
        :param sample_info: dict
        :param experimental_data: dict
        :param processed_data: dict
        :param root_dict: dict
        """
        self.sample_info= sample_info
        self.experimental_data= experimental_data
        self.processed_data= processed_data
        self.root_dict= root_dict
        self.btn_experiment_status_2.setText("Safe")
        self.btn_experiment_status_2.setStyleSheet("background-color:green")
        print("#"*100)
        print("#" * 100)
    def new_measure_of_data(self):
        """
        Manage another measurement of the sample at a given wavelength
        """
        for selected_sample in self.ls_list_of_all_samples.selectedItems():
            sample = selected_sample.data(0x0100)
        for selected_wavelength in self.ls_done_wavelengths.selectedItems():
            wavelength = selected_wavelength.data(0x0100)

        logging.info(f"Measuring again specific data of {selected_sample} at {selected_wavelength}.")

        list_of_measure_to_do = []
        for selected_data in self.ls_measure_data.selectedItems():
            a= self.ls_measure_data.row(selected_data)
            list_of_measure_to_do.append(int(a))

        logging.info(f"List of the new measure: {list_of_measure_to_do}.")

        self.measurement(sample=sample, wavelength=wavelength, new_measure=True, data_for_new_measure=list_of_measure_to_do)

        self.tpef_calculation_for_single_measure(sample=sample, wavelength=wavelength, simulation=False)
        self.btn_new_measure_data.setEnabled(True)
        self.btn_remove_data.setEnabled(True)
    def prepare_measurements(self):
        """
        Manage the samples before measurement
        :return: True
        """
        sleep(0.1)
        """Disable all widgets to avoid user modification which may cause troubles"""
        self.disable_parameters_widgets()
        """Set the order of samples for the measures"""
        self.order_of_measure.clear()
        self.order_of_measure.append(self.fluorescence_ref_dye)
        self.counter_results[self.fluorescence_ref_dye] = True
        self.sample_info[f"{self.fluorescence_ref_dye}"]["state"] = True

        """temp_dict, solvent_ref, phi_ref= rdf().reference_dye_info(self.fluorescence_ref_dye)
        wavelength_ref = []
        intensity_ref = []
        for key, value in temp_dict.items():
            wavelength_ref.append(key)
            intensity_ref.append(value)
        save().save_sample(self.fluorescence_ref_dye, wavelength_ref, intensity_ref, self.ref_concentration, solvent_ref, phi_ref)"""

        for key in self.sample_info:
            wavelength_sample = self.sample_info[key]["emission_wavelength"]
            intensity_sample = self.sample_info[key]["emission_intensity"]
            self.order_of_measure.append(key)
            save().save_sample(key, wavelength_sample, intensity_sample, self.sample_info[key]["concentration"],
                               self.sample_info[key]["solvent"], self.sample_info[key]["phi"])

        if self.first_measure == True:
            self.create_process_data_files()
            self.first_measure = False

        logging.info(f"order of measure: {self.order_of_measure}.")

        create_samples(self.order_of_measure)
        for sample in self.order_of_measure:
            save().create_file(sample)

        logging.info("initialistion ok !")

        return True
    def pop_tune_opo_wavelength(self, wavelength):
        """
        Message to user for tuning the wavelength of the OPO manually, then continue the measure when oko is pressed.
        :param wavelength: int, excitation wavelength
        :return: bool
        """

        logging.info("OPO tuning")

        dlg = QMessageBox(self)
        dlg.setWindowTitle(f"OPO wavelength tuning")
        dlg.setText(f"Laser is running ! \nPlease tune the wavelength of the OPO to {wavelength} nm, then when ready for measure press OK.")
        dlg.exec()
        if QMessageBox.StandardButton.Ok:
            return True
        elif QMessageBox.StandardButton.Cancel:
            return False
    def pop_up(self, name):
        """
        Send a message to user to change sample
        :param name: str, sample name
        """
        dlg = QMessageBox(self)
        dlg.setWindowTitle(f"Next sample is {name}")
        dlg.setText(f"Please put {name} in the sample holder !")
        dlg.exec()
        if QMessageBox.StandardButton.Ok:

            logging.info(f"User as changed sample to {name}.")

    def pop_up_end_of_measure(self):
        """
        Message send to the user to say that the measurement of the sample at this wavelength is finished.
        """
        self.plot_s2_from_signal(simulation=False)
        self.btn_run_tpa.setEnabled(True)
        dlg = QMessageBox(self)
        dlg.setWindowTitle(f"The measurement is done.")
        dlg.setText(f"The measurement is done.")
        dlg.exec()
        if QMessageBox.StandardButton.Ok:
            return
    def run_tpa_measurements(self):
        """
        Check that everything is OK before measurement to avoid crashs.
        """
        self.change_progress_bar_measurment_status(0)
        for key, value in self.checklist.items():

            logging.info(f"self.checklist_key:{key}, value:{value}.")

            if value == False:
                self.ls_terminal_parameters.clear()
                self.ls_terminal_parameters.addItem(f"The {key} is not done")

                logging.warning(f"The {key} is not done.")

                return
            else:
                corr_dye_area = rdf()._Read_Reference_Fluo(self.fluorescence_ref_dye, emission_filter=self.current_emission_filter)

                logging.info(f"Correction of the reference dye finish, value= {corr_dye_area}.")

                spectrum, solvent, phi = rdf().reference_dye_info(self.fluorescence_ref_dye)
                self.save_sample_informations(phi=phi, solvent= solvent)
                emission_wavelength = []
                emission_intensity = []
                for key, value in spectrum.items():
                    emission_wavelength.append(key)
                    emission_intensity.append(value)
                self.sample_info[f"{self.fluorescence_ref_dye}"] = {"emission": None,
                                                               "concentration": self.ref_concentration,
                                                               "solvent": solvent, "phi": phi,
                                                               "emission_wavelength": emission_wavelength,
                                                               "emission_intensity": emission_intensity,
                                                               "Area correction factor": corr_dye_area,
                                                               "CCD reference": self.ccd_correction_dye,
                                                               "state": True}
            """for key in self.sample_info:
                print(f"Informations of {key} are ", self.sample_info[key])
            print("Ready for sample measurement !")
            print("Starting initialisation !")"""
            preparation = self.prepare_measurements()

            logging.info(f"{preparation}")

        if preparation == False:

            logging.info("Not all parameters are given, please fill them, then restart")

            return
            """If user don't choose any, the first sample & wavelength will be measured"""
        else:
            sample= self.select_sample()
            wavelength= self.select_wavelength()
            """Here the user choose which sample & wavelength to measure"""

            logging.info("data refreshed")
            logging.info(f"Ready for measurement of {sample} at {wavelength} nm.")

            if wavelength > 1080:
                self.opo_measure_state=True
                fs_laser().setWavelengthBlocking(800)
                fs_laser().openShutterBlocking()

                logging.info(f"Laser tuned to {wavelength}.")
                logging.info("Laser shutter is open.")

                self.pop_tune_opo_wavelength(wavelength=wavelength)
                fs_laser().closeShutterBlocking()
                self.opo_wavelength=wavelength
                wavelength = 800
                print(self.opo_wavelength)
                self.measurement(sample=f"{sample}", wavelength=wavelength, new_measure=False,
                                 data_for_new_measure=None, opo=True)
            else:
                self.opo_measure_state=False
                self.opo_wavelength=None
                self.measurement(sample= f"{sample}", wavelength= wavelength, new_measure= False, data_for_new_measure=None)
    def save_fluo_vs_power(self, i, temporary_wavelength_list, temporary_intensity_list, sample, wavelength):
        """
        Save fluorescence vs power file as png.
        :param i: int, measure n°
        :param temporary_wavelength_list: array, wavelength
        :param temporary_intensity_list: array, intensity
        :param sample: str, sample name
        :param wavelength: int, excitation wavelength
        """
        save().figure_intensity_vs_power(number_of_measure=i, x=temporary_wavelength_list,
                                         y=temporary_intensity_list, sample=sample, wavelength=wavelength)
    def save_quad_from_signal(self, square_power, fluo, fit, sample, wavelength):
        """
        Save quadraticity check file as png.
        :param square_power: array, P²
        :param fluo: array, F
        :param fit: array, slope of log(F)/log(P²)
        :param sample: str, sample name
        :param wavelength: int, excitation wavelength
        """
        save().figure_quadraticity(x1=square_power, y1=fluo, x2=square_power, y2=fit,
                                   sample=sample, wavelength=wavelength)
    def save_plot_from_signal(self, x, y, wavelength, sample):
        """
        Save dark noise file as png.
        :param x:array, wavelength
        :param y: array, intensity
        :param wavelength: int, excitation wavelength
        :param sample: str, sample name
        """
        save().figure(x=x, y=y, name=f"Dark at {wavelength}", sample=sample)
    def save_sample_informations(self, solvent, phi):
        """
        This function save data from user for post treatment.
        :param solvent: str
        :param phi: float
        """
        number_of_measures= self.power_pitch+2
        data_to_save= {'CCD':self.ccd_correction_dye,
                       'TPEF':{
                           "name": self.fluorescence_ref_dye,
                           "concentration": self.ref_concentration,
                           "phi":phi,
                           "solvent": solvent
                               },
                       'i':number_of_measures,
                       'Emission_filter':self.current_emission_filter,
                       'Integration time':self.integration_time_value
                       }
        for key in self.sample_info:
            data_to_save.update({f"{key}":{
                                            "phi":self.sample_info[key]["phi"],
                                            "concentration":self.sample_info[key]["concentration"],
                                            "solvent": self.sample_info[key]["solvent"]
                                            }
                                })
        save().create_samples_informations(data_to_save=data_to_save)

    """
    Measurement _  Processing data
    """
    def clean_labels(self):
        """
        Clear ls_measure_data widget
        """
        self.ls_measure_data.clear()
    def clear_plot_from_signal(self):
        """
        Clear PlotWidgets S2Spectrum, FluoIntensitySpectrum, FPSpectrum_2.
        """

        logging.info("Clear plotwidgets.")

        self.S2Spectrum.clear()
        self.FluoIntensitySpectrum.clear()
        self.FPSpectrum_2.clear()
        self.S2Spectrum.addLegend()
        self.FluoIntensitySpectrum.addLegend()
    def get_files(self):
        """
        Read sample informations for the selected sample and update processed data dictionary.
        """
        self.sample_info = sdp().read_samples_informations(root=self.folder_root)

        logging.info(f"self.sample_info: {self.sample_info}.")

        samples= get_samples()

        logging.info(f"sample: {samples}.")

        self.list_of_sample_data ={}
        self.list_of_wavelength=[]
        self.list_of_sample_data[samples[0]] = sdp().read_process(samples[0], simulation=False, root=self.folder_root)
        for wave in self.list_of_sample_data[samples[0]]:
            self.list_of_wavelength.append(int(wave))

            logging.info(f"self.list_of_waveelngth: {self.list_of_wavelength}, type: {type(self.list_of_wavelength)}.")

        self.refresh_wavelength(self.list_of_wavelength)
        n, data = sdp().extract_data_from_experimental(sample=samples[0], wavelength=self.list_of_wavelength[0], root=self.folder_root)
        i=1
        for number_of_measure in n:
            i+=1

            logging.info(f"measure n°: {number_of_measure}.")
        logging.info(f"Total measures= {i}.")

        self.power_pitch = i
        self.fluorescence_ref_dye = self.sample_info["TPEF"]["name"]

        logging.info(f"self.fluorescence_dye: {type(self.fluorescence_ref_dye)}.")

        for sample in samples:

            logging.info(f"Sample: {sample}.")

            self.processed_data.update({f"{sample}":sdp().read_process(sample= sample, simulation=False, root=self.folder_root)})
            sdp().copy_processed_data_for_simulation(sample, root=self.folder_root)
    def graph_data(self, x, y, view, order_process, x_log, y_log, x_log2, y_log2):
        """
        Plot data on view.
        :param x: array ou list
        :param y: array ou list
        :param view: str
        :param order_process: float
        :param x_log: array
        :param y_log: array
        :param x_log2: array
        :param y_log2: array
        """
        x = np.array(x)
        y = np.array(y)
        if view == "fluo":
            self.FluoIntensitySpectrum.plot(x, y, show= True)
        elif view == "S2F":
            self.S2Spectrum.clear()
            x, y= self.sort_wavelength(x, y)
            self.S2Spectrum.plot(x, y, show = True)
        elif view == "quad":
            self.FPSpectrum_2.clear()
            self.FPSpectrum_2.addLegend()
            self.FPSpectrum_2.setLabel('left', text= "log(Fluorescence intensity) / cps")
            self.FPSpectrum_2.setLabel('bottom', text="log(P)")
            self.FPSpectrum_2.plot(x_log, y_log, pen= "k", symbol= "o", symbolPen="g",symbolBrush=0.5, show= True)
            self.FPSpectrum_2.plot(x_log2, y_log2, show=True)
            self.FPSpectrum_2.setLabel('top', text=f"Order of the absorption = {order_process}")
    def plot_from_signal(self, x, y, square_power_dark_corr, full_corrected_area, i, power):
        """
        Plot all fluorescence spectrum vs laser power
        :param x: array ou list
        :param y: array ou list
        :param square_power_dark_corr: array
        :param full_corrected_area: array
        :param i: int
        :param power: float
        :return:
        """

        logging.info(f"square_power_dark_corr: {square_power_dark_corr}, full_corrected_area: {full_corrected_area}.")

        dictionary_color = {"0": "w", "1": "b", "2": "c", "3": "g", "4": "y", "5": "m", "6": "r", "7": "w",
                            "8": "b", "9": "c", "10":"g"}

        try:
            power= int(power)
        except ValueError:

            logging.warning('Cannot convert power as a integer.')

        x = np.array(x)
        y = np.array(y)
        self.FluoIntensitySpectrum.addLegend()
        self.FluoIntensitySpectrum.setLabel('left', text="Fluorescence intensity / cps")
        self.FluoIntensitySpectrum.setLabel('bottom', text="Wavelength / nm")
        if i > 10:

            logging.info(f"i: {i}.")

            data=str(i)
            data_temp=data[1:]
            i=int(data_temp)

            logging.info(f"i: {i}.")

        self.FluoIntensitySpectrum.plot(x, y, name= f"{power} \n", pen= dictionary_color[f"{i}"], show=True)
    def plot_s2_from_signal(self, simulation):
        """
        Plot S2 values on PlotWidget S2Spectrum
        :param simulation: bool
        """
        try :
            self.S2Spectrum.clear()
            self.show_data_of_sample(simulation=simulation)
            sample = self.select_sample()
            data = sdp().read_process(sample=sample, simulation=simulation, root=self.folder_root)
            sigma_2 = []
            tpef_wavelength = []
            for wavelength in data:
                sigma_2.append(int(data[f"{wavelength}"]["S2"]))

                logging.info(f"Type sigma2: {type(sigma_2)}.")

                tpef_wavelength.append(int(wavelength))

                logging.info(f"Type wavelength: {type(wavelength)}.")

            x = np.array(tpef_wavelength)
            y = np.array(sigma_2)
            x, y= self.sort_wavelength(x, y)
            logging.info(f"Wavelength: {x}.")
            logging.info(f"S2: {y}.")

            self.S2Spectrum.setLabel('left', text="Sigma 2 / GM")
            self.S2Spectrum.setLabel('bottom', text="Wavelength / nm")
            self.S2Spectrum.plot(x, y, pen="g", symbol="o", symbolPen="m", symbolBrush=0.5, name=f"{sample}", show=True)
            sdp().figure(x=x, y=y, name=f"TPEF spectra of {sample}", sample=sample, root=self.folder_root)
        except:

            logging.warning("error in plot graph")

    def refresh_data(self):
        """
        Refresh send data
        """
        self.ls_list_of_all_samples.clear()
        samples = get_samples()
        for sample in samples:
            lw_item = QListWidgetItem(sample.sample_name)
            lw_item.setData(0x0100, sample)  # 0x0100 => QtCore.Qt.UserRole
            self.ls_list_of_all_samples.addItem(lw_item)

            logging.info(f"sample: {sample}.")

    def refresh_wavelength(self, list_of_wavelength):
        """
        Refresh wavelength values in ls_done_wavelength
        :param list_of_wavelength:list
        """
        self.ls_done_wavelengths.clear()
        for wavelength in list_of_wavelength:
            wavelength= int(wavelength)
            lw_item = QListWidgetItem(f"{wavelength}")
            lw_item.setData(0x0100, wavelength)  # 0x0100 => QtCore.Qt.UserRole
            self.ls_done_wavelengths.addItem(lw_item)

            logging.info("User refresh wavelength list.")
            logging.info(f"Wavelength= {wavelength}.")

    def read_data_from_folder(self):
        """
        Read data of samples from files samples.json and sample_data.json
        """

        logging.info("Reading data of samples.")

        temp = cfg.get_path()
        self.folder_root = Path(temp)
        self.sample_root = self.folder_root / "Metadata" / "samples.json"
        self.sample_data_root = self.folder_root / "Metadata" / "sample_data.json"

        logging.info(f"Folder_root: {self.folder_root}.")
        logging.info(f"Sample_data_root: {self.sample_data_root}.")

        self.refresh_data()
        self.get_files()
    def remove_data(self):
        """
        Remove selected data from the experimental data files
        """

        logging.info("User removed data.")

        list_of_measure_to_remove = []
        for selected_data in self.ls_measure_data.selectedItems():
            a= self.ls_measure_data.row(selected_data)

            logging.info(f"Data to remove: {a}.")

            list_of_measure_to_remove.append(int(a))

        logging.info(f"List_of_measure_to_remove: {list_of_measure_to_remove}.")

        for selected_sample in self.ls_list_of_all_samples.selectedItems():
            sample = selected_sample.data(0x0100)
        for selected_wavelength in self.ls_done_wavelengths.selectedItems():
            wavelength = selected_wavelength.data(0x0100)
        for i in list_of_measure_to_remove:

            logging.info(f"Delete: #{i}, sample: {sample}, wavelength: {wavelength}.")

            sdp().remove_data(sample=sample, wavelength=wavelength, measure=i, simulation=False,root=self.folder_root)
        self.tpef_calculation_for_single_measure(sample=sample, wavelength=wavelength, simulation= False)
        self.plot_s2_from_signal(simulation=False)

        logging.info("Data removed => ok !")

    def save_fluo_vs_power(self, i, temporary_wavelength_list, temporary_intensity_list, sample, wavelength):
        """
        Save fluorescence vs laser power as png file
        :param i: int
        :param temporary_wavelength_list: list
        :param temporary_intensity_list: list
        :param sample: str
        :param wavelength: int
        """
        sdp().figure_intensity_vs_power(number_of_measure=i, x=temporary_wavelength_list,
                                         y=temporary_intensity_list, sample=sample, wavelength=wavelength, root=self.folder_root)

        logging.info("save_fluo_vs_power => ok !")

    def save_plot_from_signal(self, x, y, wavelength, sample):
        """
        Save dark noise as png file
        :param x: array
        :param y: array
        :param wavelength: int
        :param sample: str
        """
        sdp().figure(x=x, y=y, name=f"Dark at {wavelength}", sample=sample, root=self.folder_root)

        logging.info("save_plot_from_signal => ok !")

    def save_quad_from_signal(self, square_power, fluo, fit, sample, wavelength):
        """
        save quadraticity check as png file.
        :param square_power: array
        :param fluo: array
        :param fit: array
        :param sample: str
        :param wavelength: int
        """
        sdp().figure_quadraticity(x1=square_power, y1=fluo, x2=square_power, y2=fit,
                                   sample=sample, wavelength=wavelength, root=self.folder_root)

        logging.info("save_plot_from_signal => ok !")

    def select_sample(self):
        """
        Return selected sample as str
        :return: str
        """
        for selected_sample in self.ls_list_of_all_samples.selectedItems():
            sample = selected_sample.data(0x0100)

            logging.info(f"The selected sample is : {sample}.")

            return sample
    def select_wavelength(self):
        """
        Return selected wavelength as str
        :return: str
        """
        for selected_wavelength in self.ls_done_wavelengths.selectedItems():
            wavelength = selected_wavelength.data(0x0100)

            logging.info(f"The selected wavelength is: {wavelength}.")

        return wavelength
    def show_data_of_sample(self, simulation):
        """
        Show experimental data on ls_measure_data widget of selected sample and wavelength.
        :param simulation: bool
        """

        logging.info("Starting show_data_of_sample.")

        self.clean_labels()
        number_of_measure = self.power_pitch+2
        for selected_sample in self.ls_list_of_all_samples.selectedItems():
            sample = selected_sample.data(0x0100)

            logging.info(f"Selected sample: {sample}.")

        for selected_wavelength in self.ls_done_wavelengths.selectedItems():
            wavelength = selected_wavelength.data(0x0100)

            logging.info(f"Selected wavelength: {wavelength}.")

            if simulation == False:

                logging.info("Simulation is false.")

                file_root = self.folder_root / f"{sample}" / f"experimental_data {sample}_{wavelength}.json"
                list_of_fl, list_of_pp = sdp()._get_quick_f_and_pp(sample=sample, wavelength=wavelength,
                                                                    root=self.folder_root, number_of_measure=number_of_measure)

                logging.info("_get_quick_f_and_pp => ok !")

                data = sdp().read_process(sample=sample, simulation=simulation, root=self.folder_root)

                logging.info("read_process => ok !")

                slope= data[f"{wavelength}"]["quadraticity"]['slope']
                coeff=data[f"{wavelength}"]["quadraticity"]['coeff']
                log_square_fluo=data[f"{wavelength}"]["quadraticity"]["log(F)"]
                log_square_power=data[f"{wavelength}"]["quadraticity"]['log(P²)']
                fit=data[f"{wavelength}"]["quadraticity"]['fit']
                full_corrected_area=data[f"{wavelength}"]['full_corrected_area']
                sigma2_phi=data[f"{wavelength}"]['S2F']
                sigma2=data[f"{wavelength}"]['S2']
                data.update(
                        {
                            f"{wavelength}":
                                {
                                    "quadraticity":
                                        {
                                            'slope': slope,
                                            'coeff': coeff,
                                            'log(F)': log_square_fluo,
                                            'log(P²)': log_square_power,
                                            'fit': fit,
                                            'F': list_of_fl,
                                            'PP': list_of_pp
                                        },
                                    'full_corrected_area': full_corrected_area,
                                    "S2F": sigma2_phi,
                                    "S2": sigma2
                                }
                        }
                )

                logging.info("Data dictionnary updated => ok !")

                sdp().save_process(dictionary=data, sample=sample, state=False,
                                               simulation=simulation, root=self.folder_root)

                logging.info("save_processs => ok !")

            elif simulation == True:

                logging.info("Simulation is True.")

                file_root = self.folder_root / f"{sample}" / f"experimental_data {sample}_{wavelength}_temp.json"
            data = sdp().read_json_data(file_root)

            logging.info(f"File root : {file_root}.")

            for i in range(number_of_measure):
                try:

                    logging.info(f"Data dictionnary keys: {data.keys()}")

                    fluorescence_intensity = data[f"{wavelength}"][f"{i}"]['fully corrected area']
                    power = data[f"{wavelength}"][f"{i}"]['raw power - dark']
                    square_power = power*power
                    quadraticity = (float(fluorescence_intensity)/float(square_power))
                    square_power = round(square_power, 8)

                    logging.info(f"Power: {power}, quadraticity: {quadraticity}, square power: {square_power}.")

                    self.ls_measure_data.addItem(f"{i} \t {fluorescence_intensity:.2E} \t {square_power} \t {quadraticity:.2E}")
                except KeyError:

                    logging.error("key error _show data sample !")

        s2f, x_wavelength, log_f, log_pp, list_fl, list_pp = sdp().extract_process_data(sample=sample,
                                                                       wavelength=wavelength, simulation=simulation, root=self.folder_root)

        try:
            x2, y2, value, order_process, x_log, y_log, x_log2, y_log2= sdp().fit_curve(x= np.array(list_pp), y=np.array(list_fl))

            self.graph_data(x=list_pp, y=list_fl, order_process= order_process, view="quad", x_log=x_log, y_log=y_log, x_log2=x_log2, y_log2=y_log2)

            logging.info("Fitting curve => ok !")

        except:

            logging.warning("Impossible to fit")

        self.FluoIntensitySpectrum.clear()
        i=0
        data_bis= sdp().extract_data_from_experimental_bis(sample= sample, wavelength= wavelength, root= self.folder_root)
        for key in data_bis:
                for power in data_bis[f"{key}"]:
                    if power != "dark":
                        i += 1
                        try:

                            logging.info(f"Power= {power}.")

                            x, y = sdp()._get_intensity_vs_wavelength_data(data_bis[f"{key}"][f"{power}"])

                            logging.info(f"Length x: {len(x)}")
                            logging.info(f"Length y: {len(y)}")

                            self.plot_from_signal(x= x, y=y, square_power_dark_corr="", full_corrected_area="", i= i, power= i)
                        except:
                            print("Error for plot")

    def show_data_to_user(self, data):
        """
        Communication between class Measure and user
        :param data: str
        """
        self.communication_with_user(data, clear=False, measurement=True)
    def simulation(self):
        """
        Manage the simulation of fitting log(F)/log(P²)
        """

        logging.info("Simulation running")

        list_of_measure_to_remove = []
        for selected_data in self.ls_measure_data.selectedItems():
            a = self.ls_measure_data.row(selected_data)

            logging.info(f"Selected data is: {selected_data}.")

            list_of_measure_to_remove.append(a)
        for selected_sample in self.ls_list_of_all_samples.selectedItems():
            sample = selected_sample.data(0x0100)

            logging.info(f"Selected sample is: {selected_sample}.")

        for selected_wavelength in self.ls_done_wavelengths.selectedItems():
            wavelength = selected_wavelength.data(0x0100)

            logging.info(f"Selected wavelength is: {selected_wavelength}.")

        sdp().copy_data_for_simulation(sample=sample, wavelength= wavelength, root=self.folder_root)
        sdp().copy_processed_data_for_simulation(sample=sample, root=self.folder_root)
        for i in list_of_measure_to_remove:

            logging.info(f"List of measure to remove: {list_of_measure_to_remove}.")
            logging.info(f"simulate without : #{i}, sample: {sample}, wavelength: {wavelength}.")

            sdp().remove_data(sample=sample, wavelength=wavelength, measure=i, simulation=True, root=self.folder_root)

        self.tpef_calculation_for_single_measure(sample=sample, wavelength=wavelength, simulation=True)
        self.plot_s2_from_signal(simulation=True)

        logging.info("Simulation => ok !")

    def sort_wavelength(self, x, y):
        """
        Sort the data for to get a cleaner plot
        :param x: array
        :param y: array
        :return: array, array
        """
        size_of_data = x.size
        dict = {}
        for n in range(size_of_data):
            dict[x[n]]=y[n]
        print(dict)
        x_sort=sorted(x)
        y_temp = []
        for n in range(size_of_data):
            a= dict[x_sort[n]]
            y_temp.append(a)
        print(y_temp)
        print(x_sort)
        y=np.array(y_temp)
        return x_sort, y





    def tpef_calculation_for_single_measure(self, sample, wavelength, simulation):
        """
        Caclulate the Sigma 2 value of the sample at selected wavelength
        :param sample: str
        :param wavelength: int
        :param simulation:  bool
        """
        list_of_fully_corrected_area= []
        list_of_square_dark_corr_power = []
        data = sdp()._get_data_for_new_mesure(sample= sample, wavelength= wavelength, simulation=simulation, root=self.folder_root)
        for key_wavelength in data:

            logging.info("Starting TPEF calulation.")
            logging.info(f"Type data : {type(data)}.")

            for key_measure in data[f"{key_wavelength}"]:
                if key_measure != "dark":

                    logging.info(f"key_measure: {key_measure}.")
                    logging.info(f"type(key_measure): {type(key_measure)}.")

                    power = data[f"{key_wavelength}"][f"{key_measure}"]["raw power - dark"]
                    square_power = power*power

                    logging.info(f"Square power = {square_power}.")

                    list_of_square_dark_corr_power.append(square_power)
                    list_of_fully_corrected_area.append(data[f"{key_wavelength}"][f"{key_measure}"]["fully corrected area"])

        area, correlation_coeff, slope, coeff, log_square_fluo, log_square_power, fit = rdf()._Quad_Log(list_of_fully_corrected_area,
                                                                               list_of_square_dark_corr_power)


        logging.info(f"self.fluorescence_ref_dye: {self.fluorescence_ref_dye}, type(self.fluorescence_ref_dye): {type(self.fluorescence_ref_dye)}.")
        logging.info(f"sample: {sample}, type(sample): {type(sample)}.")

        sigma2_phi, sigma2 = rdf_2()._Calcul_TPA(sample_info=self.sample_info,
                                               processed_data=self.processed_data,
                                               sample_name=str(sample),
                                               reference_name=self.fluorescence_ref_dye,
                                               wavelength=int(wavelength),
                                               simulation=simulation,
                                               root=self.folder_root)
        logging.info("Update dictionnary")

        self.update_dictionary_processed_data(sample=sample, wavelength=wavelength, slope=slope,
                                              coeff=coeff, log_square_fluo=log_square_fluo,
                                              log_square_power=log_square_power, fit=fit,
                                              full_corrected_area=area,
                                              sigma2_phi=sigma2_phi, sigma2=sigma2, list_of_fl=list_of_fully_corrected_area, list_of_pp=list_of_square_dark_corr_power)

        #logging.info(f"self.processed_data: {self.processed_data}.")

        sleep(0.1)

        file_process = sdp().save_process(dictionary=self.processed_data[sample], sample=sample, state=False, simulation=simulation, root=self.folder_root)
        self.root_dict.update({sample: file_process})
        if simulation == False:

            logging.info("simulation is false! ")

            file_exp = sdp()._get_file_for_new_mesure(sample= sample, wavelength= wavelength, root= self.folder_root)
            sdp().save_experimental_to_excel(file=file_exp, sample=sample, root= self.folder_root)
            sdp().save_process_data_to_excel(file=file_process, sample=sample, root= self.folder_root)
            self.root_dict.update({sample: file_process})

            file = self.root_dict[sample]
            data_for_sigma_2 = sdp().read_json_data(file)

            logging.info(f"sample: {sample}.")
            logging.info(f"data_for_sigma_2.keys(): {data_for_sigma_2.keys()}.")

            x = []
            y = []
            for wavelength in self.wavelength_list:

                logging.info(f"data_for_sigma_2.keys(): [data_for_sigma_2.keys()].")
                logging.info(f"data_for_sigma_2.values: {data_for_sigma_2.values}.")

                x.append(wavelength)
                y.append(data_for_sigma_2[f"{wavelength}"]["S2F"])
            sdp().figure(x=x, y=y, name=f"TPEF spectra of {sample}", sample=sample, root= self.folder_root)
    def update_dictionary_processed_data(self, sample, wavelength, slope, coeff, log_square_fluo, log_square_power, fit,
                                         full_corrected_area, sigma2_phi, sigma2, list_of_fl, list_of_pp):
        """
        Update processed data dictionary
        :param sample: str
        :param wavelength: int
        :param slope: flaot
        :param coeff: float
        :param log_square_fluo: array
        :param log_square_power: array
        :param fit: float
        :param full_corrected_area:float
        :param sigma2_phi: float
        :param sigma2: float
        :param list_of_fl: array or list
        :param list_of_pp: array or list
        """
        self.processed_data.update({
            sample:
                {
                    f"{wavelength}":
                        {
                            "quadraticity":
                                {
                                    'slope': slope,
                                    'coeff': coeff,
                                    'log(F)': log_square_fluo,
                                    'log(P²)': log_square_power,
                                    'fit': fit,
                                    'F': list_of_fl,
                                    'PP': list_of_pp
                                },
                            'full_corrected_area': full_corrected_area,
                            "S2F": sigma2_phi,
                            "S2": sigma2
                        }
                }
        })

        logging.info("Processed data dictionary Updated.")



    """
    Advanced Calibration
    """
    """def reset_power_angle_conversion_file(self):
        tuner.RotationMount().measure_power_angle_conversion()
        print("Reset of the conversion file finsihed")"""

    def select_filter_E650sp2p(self):
        """
        Select emission filter
        :return: True
        """
        rdf().select_emission_filter(emission_filter="E650SP-2P")
        self.current_emission_filter = "E650SP-2P"

        logging.info(f"User select {self.current_emission_filter}.")

        return True

    def select_filter_E750sp2p(self):
        """
        Select emission filter
        :return: True
        """
        rdf().select_emission_filter(emission_filter="E750SP-2P")
        self.current_emission_filter = "E750SP-2P"

        logging.info(f"User select {self.current_emission_filter}.")

        return True

    def select_filter_E650(self):
        """
        Select emission filter
        :return: True
        """
        rdf().select_emission_filter(emission_filter="E650")
        self.current_emission_filter = "E650"

        logging.info(f"User select {self.current_emission_filter}.")

        return True

    def select_filter_FES700(self):
        """
        Select emission filter
        :return: True
        """
        rdf().select_emission_filter(emission_filter="FES_700")
        self.current_emission_filter = "FES_700"

        logging.info(f"User select {self.current_emission_filter}.")

        return True

    def select_filter_FES750(self):
        """
        Select emission filter
        :return: True
        """
        rdf().select_emission_filter(emission_filter="FES_750")
        self.current_emission_filter = "FES_750"

        logging.info(f"User select {self.current_emission_filter}.")

        return True

    def select_filter_FES800(self):
        """
        Select emission filter
        :return: True
        """
        rdf().select_emission_filter(emission_filter="FES_800")
        self.current_emission_filter = "FES_800"

        logging.info(f"User select {self.current_emission_filter}.")

        return True

    def select_filter_FF01650(self):
        """
        Select emission filter
        :return: True
        """
        rdf().select_emission_filter(emission_filter="FF01-650SP25")
        self.current_emission_filter = "FF01-650SP25"

        logging.info(f"User select {self.current_emission_filter}.")

        return True


if __name__ == "__main__":
    app = QApplication(sys.argv)
    app.setStyleSheet('''
        QWidget {
        font-size: 30 px;
        }
    ''')

    myApp = MyApp()
    myApp.show()
    sys.exit(app.exec())

    try:
        sys.exit(app.exec())
    except SystemExit:
        print('Closing Windows')