"""
Created on 5/1/2023
@author: Jonathan DANIEL, jonathan.daniel@u-bordeaux.fr
"""
import copy
from pathlib import Path
import os.path
import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
from scipy.stats import linregress
import csv
import time
from Ressources_scripts import References_data
from Ressources_scripts.configmanagement import save
import logging

logging.basicConfig(level=logging.INFO,
                    filename="PyTPEF.log",
                    filemode="a",
                    format='%(asctime)s - %(levelname)s - %(message)s')

CUR_DIR = os.path.dirname(__file__)
ROOT = Path(f"{CUR_DIR}").parent
ROOT = ROOT / "Ressources" / "save_directory.json"

EMISSION_FILTER={"FES_750":(542.0,684.0), "FES_800":(397.0,799.0), "FES_700":(433.0,695.0), "E750SP-2P":(509.9, 750.1), "E650SP-2P":(431.0,760.0), "E650":(480.5, 625.0), "FF01-650SP25":(429.9,630.1)}

class Read_Data_File():

    def __init__(self):
        """Dictionaries"""
        self.refractive_index_dic = {"Water": 1.333, "Toluene": 1.4961, "THF": 1.405, "Methanol": 1.3288, "DMSO": 1.477}
        self.spectrum_JOD_water = References_data.REF_JOD_WATER_FLUO
        self.spectrum_DMANs_toluen = References_data.REF_DMANS_TOLUEN_FLUO
        self.spectrum_fluo_fluorescein_water = References_data.REF_FLUORESCEIN_WATER_FLUO
        self.spectrum_fluo_NR_DMSO = References_data.REF_NILERED_DMSO_FLUO
        self.spectrum_S2F_NR_DMSO = References_data.REF_NILERED_DMSO_S2F
        self.spectrum_S2F_fluorescein_water =References_data.REF_FLUORESCEIN_WATER_S2F

        """Lists"""
        self.dye_fluo = []
        self.dye_limited_fluo = []
        self.dye_wavelength = []
        self.dye_limited_wavelength = []


    def select_emission_filter(self, emission_filter):
        """
        This function update the corrections data depending on the selected emission filter
        :param emission_filter: str
        :return: bool
        """
        data=EMISSION_FILTER[f"{emission_filter}"]
        self.minimum_wavelength_for_correction=data[0]
        self.maximum_wavelength_for_correction=data[1]

        logging.info(f"Tpa_calc:Select_emission_filter: {emission_filter}.")
        logging.info(f"Tpa_calc:Min_emission correction= {self.minimum_wavelength_for_correction}, Max_emission correction= {self.maximum_wavelength_for_correction}.")

        return True

    def _Calcul_TPA(self, sample_info, processed_data, sample_name, reference_name, wavelength):
        """
        This function calculates the two-photon cross section sigma 2 of sample at a given wavelength.
        :param sample_info: dict
        :param processed_data: dict
        :param sample_name: str
        :param reference_name: str
        :param wavelength: int
        :return: float, float
        """

        logging.info(f"Tpa_calc:Starting _Calcul_TPA")
        logging.info(f"Tpa_calc:sample_info: {sample_info}.")
        logging.info(f"Tpa_calc:type(processed_data): {type(processed_data)}.")
        logging.info(f"Tpa_calc:sample_name: {sample_name}.")
        logging.info(f"Tpa_calc:reference_name: {reference_name}.")
        logging.info(f"Tpa_calc:wavelength: {wavelength}.")

        if sample_name == reference_name:
            ref_data= processed_data[sample_name]
        else:
            try:
                ref_data= save().read_process(sample= reference_name)

                logging.info(f"Tpa_calc:Ref_data: {ref_data}.")

            except FileNotFoundError:

                logging.error("Tpa_calc:File not Found. Please measure the reference before the sample !")

                return -1, -1

        try:

            logging.info(f"Tpa_calc:sample_name: {sample_name}.")
            logging.info(f"Tpa_calc:processed_data.keys(): {processed_data.keys()}.")
            #temp_var=processed_data[f"{sample_name}"]
            #logging.info(f"Tpa_calc:processed_data[sample_name]: {temp_var}.")
            #temp_var=processed_data[f"{reference_name}"]
            #logging.info(f"Tpa_calc:processed_data[reference_name]: {temp_var}.")

        except KeyError:

            logging.error("Tpa_calc:Key error")

        phi = float(sample_info[f"{sample_name}"]["phi"])
        sample_concentration = sample_info[f"{sample_name}"]['concentration']
        ref_concentration = sample_info[f"{reference_name}"]['concentration']
        solvent = sample_info[f"{sample_name}"]["solvent"]
        ref_solvent= sample_info[f"{reference_name}"]["solvent"]
        sample_refractive_index = self.refractive_index_dic[solvent]
        reference_refractive_index = self.refractive_index_dic[ref_solvent]
        if reference_name == "fluorescein":
            reference_sigma2 = self.spectrum_S2F_fluorescein_water
        elif reference_name == "NR":
            reference_sigma2 = self.spectrum_S2F_NR_DMSO
        try:
            area = (processed_data[f"{sample_name}"][f"{wavelength}"]['full_corrected_area']) / (ref_data[f"{wavelength}"]['full_corrected_area'])

            logging.info(f"Tpa_calc:Reference for Sigma 2 calculation: {reference_name}.")
            logging.info(f"Tpa_calc:Area: {area}.")

        except KeyError:

            logging.error("Tpa_calc:Key Error. Please measure the reference before the sample !")

            return -1, -1
        concentration = float(ref_concentration) / float(sample_concentration)
        refractive_index = sample_refractive_index / reference_refractive_index
        calcul_1 = area * concentration * refractive_index

        logging.info(f"Tpa_calc:reference_sigma2.keys(): {reference_sigma2.keys()}.")

        time.sleep(0.1)
        calcul_tpa = calcul_1 * (reference_sigma2[f"{wavelength}"])

        logging.info(f"Tpa_calc:Calculation 2PA : {calcul_tpa}, type: {type(calcul_tpa)}.")
        logging.info(f"Tpa_calc:Phi: {phi}, type: {type(phi)}.")

        calcul_s2 = calcul_tpa / phi
        sigma2_phi = round(calcul_tpa, 1)
        sigma2 = round(calcul_s2, 1)

        logging.info("Tpa_calc:_Calcul_TPA => ok !")

        return sigma2_phi, sigma2

    def reference_dye_info(self, ref_dye):
        """
        This function returns informations about the reference dye.
        :param ref_dye: str
        :return: dict, str, float
        """
        if ref_dye == "fluorescein":
            solvent = "Water"
            phi = 0.9

            logging.info(f"Tpa_calc:Ref_dye: {ref_dye}, solvent: water, Phi: 0.9")

            return self.spectrum_fluo_fluorescein_water, solvent, phi
        elif ref_dye == "NR":
            solvent = "DMSO"
            phi = 0.79

            logging.info(f"Tpa_calc:Ref_dye: {ref_dye}, solvent: DMSO, Phi: 0.79")

            return self.spectrum_fluo_NR_DMSO, solvent, phi

    def _Read_Reference_Fluo(self, dye_as_reference, emission_filter):
        """
       This function reads and calculate the corrected area factor (depending on the CCD dye)
        :param dye_as_reference: str
        :param emission_filter: str
        :return: float
        """

        logging.info("Tpa_calc:Starting _Read_Reference_Fluo.")
        logging.info(f"Tpa_calc:dye_as_reference: {dye_as_reference}.")
        logging.info(f"Tpa_calc:emission_filter: {emission_filter}.")

        self.select_emission_filter(emission_filter=emission_filter)
        ref = dye_as_reference
        wavelength_min = self.minimum_wavelength_for_correction
        wavelength_max = self.maximum_wavelength_for_correction
        if ref == "fluorescein":
            spectrum = self.spectrum_fluo_fluorescein_water
        elif ref == "NR":
            spectrum = self.spectrum_fluo_NR_DMSO
        for key, value in spectrum.items():
            self.dye_wavelength.append(float(key))
            self.dye_fluo.append(float(value))
            if float(key) >= float(wavelength_min) and float(key) <= float(wavelength_max): #loop to use only wavelength inside wave_min and wave_max
                self.dye_limited_wavelength.append(float(key))
                self.dye_limited_fluo.append(float(value))
        area_full_dye = np.trapz(self.dye_fluo, self.dye_wavelength)
        area_limited_dye = np.trapz(self.dye_limited_fluo, self.dye_limited_wavelength)
        corr_dye_area = area_full_dye / area_limited_dye

        logging.info(f"Tpa_calc:area_full_dye: {area_full_dye}")
        logging.info(f"Tpa_calc:area_limited_dye: {area_limited_dye}")
        logging.info(f"Tpa_calc:corr_dye_area: {corr_dye_area}")

        self._Figure_Duo(self.dye_wavelength, self.dye_fluo, self.dye_limited_wavelength, self.dye_limited_fluo, f"{ref} limited")
        self.dye_fluo.clear()
        self.dye_limited_fluo.clear()
        self.dye_limited_wavelength.clear()
        self.dye_wavelength.clear()

        logging.info("Tpa_calc:_Read_Reference_Fluo => ok !")

        return corr_dye_area

    def _Read_Dyes_Emission_Spectrum(self, root, sample_name, emission_filter):
        """
        This function reads fluorescence intensity  files provides by user on parameters and calculate
        a correction factor of the emission area depending on the dye used for ccd correction
        :param root: path
        :param sample_name: str
        :param emission_filter: str
        :return: float
        """

        logging.info("Tpa_calc:Starting _Read_Dyes_Emission_Spectrum.")
        logging.info(f"Tpa_calc:root: {root}.")
        logging.info(f"Tpa_calc:sample_name: {sample_name}.")
        logging.info(f"Tpa_calc:emission_filter: {emission_filter}.")

        self.select_emission_filter(emission_filter=emission_filter)
        wavelength_min = self.minimum_wavelength_for_correction
        wavelength_max = self.maximum_wavelength_for_correction

        logging.info(f"Tpa_calc:wavelength_min: {wavelength_min}.")
        logging.info(f"Tpa_calc:wavelength_max: {wavelength_max}.")

        """Step 1. Read the fluorescence files which contains wavelength and corresponding fluorescence intensity
        separated by a tabulation, then put these datas in two separated lists named dye_limited_fluo or wavelength."""
        with open(root, 'r') as f:
            dye = f.readlines()  # Here I read the files of emission of the dye (column 1 : wavelength, column 2: fluo)
            for n in range(0, len(dye)):  # this is an loop with the number of wavelength
                dye_list = dye[n].split("\t")  # I copy to a temporary list and split wavelength and fluo
                self.dye_wavelength.append(float(dye_list[0]))  # first index of the list is wavelength
                self.dye_fluo.append(float(dye_list[1]))  # Second index of the list is fluorescence
                if float(dye_list[0]) >= wavelength_min and float(dye_list[0]) <= wavelength_max:
                    self.dye_limited_wavelength.append(float(dye_list[0]))
                    self.dye_limited_fluo.append(float(dye_list[1]))

        """Step 2. Determine the full area of the emission of the sample and it's area between the wavelength 
        min and max of the ccd correction dye."""
        area_full_dye = np.trapz(self.dye_fluo, self.dye_wavelength)
        area_limited_dye = np.trapz(self.dye_limited_fluo, self.dye_limited_wavelength)

        logging.info(f"Tpa_calc:area_full_dye: {area_full_dye}.")
        logging.info(f"Tpa_calc:area_limited_dye: {area_limited_dye}.")

        """Step 3. Calculate the area correction factor to takes into account the part of the fluorescence intensity
        which is not take into account by the ccd correction (out of range)."""
        corr_dye_area = area_full_dye / area_limited_dye

        logging.info(f"Tpa_calc:corr_dye_area: {corr_dye_area}.")

        """Step 4. Generate a figure with both full emission and crop emission of the sample."""
        self._Figure_Duo(self.dye_wavelength, self.dye_fluo, self.dye_limited_wavelength,
                         self.dye_limited_fluo, f"{sample_name} limited")

        """Clear all the lists for the next call."""
        self.dye_fluo.clear()
        self.dye_limited_fluo.clear()
        self.dye_limited_wavelength.clear()
        self.dye_wavelength.clear()
        """Return corr_dye_area: float
         This float will be multiplied later to ccd-corrected emission spectra of the sample."""

        logging.info("Tpa_calc:_Read_Dyes_Emission_Spectrum => ok !")

        return corr_dye_area

    def _Spectro_Correction_Factor(self, x, y_corr,emission_filter, dye_for_correction):
        """
        This function corrects the spectral response of the USB spectro using a fluorescent dye. The wavelength range
        (min & max) depends on the emission spectrum of the fluorescent dye => no correction where the dye is not
        fluorescent
        :param x: array
        :param y_corr: array
        :param emission_filter: str
        :param dye_for_correction: str
        :return: lst, float
        """

        logging.info("Tpa_calc:Starting _Spectro_Correction_Factor.")
        logging.info(f"Tpa_calc:type(x): {type(x)}.")
        logging.info(f"Tpa_calc:type(y_corr): {type(y_corr)}.")
        logging.info(f"Tpa_calc:emission_filter: {emission_filter}.")
        logging.info(f"Tpa_calc:dye_for_correction: {dye_for_correction}.")

        if dye_for_correction == "JOD":
            ref_dye = self.spectrum_JOD_water
        elif dye_for_correction == "DMANs":
            ref_dye = self.spectrum_DMANs_toluen

        self.select_emission_filter(emission_filter=emission_filter)
        wavelength_min = self.minimum_wavelength_for_correction
        wavelength_max = self.maximum_wavelength_for_correction
        wavelength_reference_tpa = []
        fluo_reference_tpa = []

        """Loop to check that the arrays will have the same size with same wavelengths, require that the emission 
        spectrum of the ccd correction dye is record with a pitch of 0.1 nm."""
        for _ in x.flat:
            _ = round(float(_), 1)
            if _ > wavelength_min and _ < wavelength_max:
                wavelength_reference_tpa.append(f"{_}")
                try:
                    if ref_dye[f"{_}"] == 0:
                        ref_dye[f"{_}"] = 1
                    fluo_reference_tpa.append(ref_dye[f"{_}"])
                except KeyError:

                    logging.error("Tpa_calc:Error in CCD corr, check filter & dye")

                    return [], True
        i = 0  # set the starting counter
        z = 0  # set the ending counter

        """This loop set the steps between wavelength_min & max."""
        for element in x.flat:
            if element <= wavelength_max:
                z += 1
            if element < wavelength_min:
                i += 1

        """Create tables of the ref dye containing only the intensity and wavelength between wavelength min & max."""
        wavelength_measured_sliced = x[i:z].copy()
        intensity_measured_sliced = y_corr[i:z].copy()
        intensity_ref = np.array(fluo_reference_tpa)
        """Check the size of the tables and normalize the correction spectrum with the minimum = 1. Later, all the
        emission spectrum are multiplied by this corrected spectrum to take into account the difference of spectral
        response of the USB spectro."""

        logging.info(f"Tpa_calc:Size of wavelength_measured_slice: {wavelength_measured_sliced.size}, First and last values: {wavelength_measured_sliced[0]} & {wavelength_measured_sliced[-1]}.")
        logging.info(f"Tpa_calc:Len of wavelength_reference_tpa: {len(wavelength_reference_tpa)}, First and last values: {wavelength_reference_tpa[0]} & {wavelength_reference_tpa[-1]}.")

        if wavelength_measured_sliced.size == len(wavelength_reference_tpa):

            logging.info(f"Tpa_calc:type(fluorescence intensity sliced): {type(intensity_measured_sliced)}.")

            corr_fact = intensity_ref / intensity_measured_sliced

            logging.info(f"Tpa_calc:type(corr_fact): {type(corr_fact)}.")
            logging.info(f"Tpa_calc:corrector factor minimum is: {corr_fact.min()}.")

            if corr_fact.min() < 0 or corr_fact.min() == 0:
                corr_fact -= (1-corr_fact.min())
                corr_fact/= (corr_fact.min())

                logging.info(f"Tpa_calc:type(corr_fact): {type(corr_fact)}.")
                logging.info(f"Tpa_calc:corrector factor minimum is: {corr_fact.min()}.")

            else:

                logging.info(f"Tpa_calc:type(corr_fact): {type(corr_fact)}.")
                logging.info(f"Tpa_calc:corrector factor minimum is: {corr_fact.min()}.")

                corr_fact /= (corr_fact.min())
        else:

            logging.error("Tpa_calc:Error in wavelength_measured_sliced.size !")
            logging.error("Tpa_calc:Corr_fact is set to 0 !")

            corr_fact = 0
        """Return corr_fact: array 
        This array correspond to a correction of spectro response between wavelength_min & max only !"""

        logging.info("Tpa_calc:_Spectro_Correction_Factor => ok !")

        return wavelength_measured_sliced, corr_fact

    def _Emission_Correction(self, corr_dye_area, table_fluo, table_wavelength, corr_fact, emission_filter):
        """Correct the emission spectrum of the sample with corr_fact (USB spectro response), then correct the emission
        area, and return a fully corrected intensity area for sigma 2 calculation)

        corr_dye_area: float
        table_fluo: array
        table_wavelength: array
        corr_fact: array

        """

        logging.info("Tpa_calc:Starting _Emission_Correction.")
        logging.info(f"Tpa_calc:corr_dye_area: {corr_dye_area}.")
        logging.info(f"Tpa_calc:table_fluo: {table_fluo}.")
        logging.info(f"Tpa_calc:table_wavelength: {table_wavelength}.")
        logging.info(f"Tpa_calc:type(corr_fact): {type(corr_fact)}.")
        logging.info(f"Tpa_calc:emission_filter: {emission_filter}.")

        self.select_emission_filter(emission_filter=emission_filter)
        wavelength_min = self.minimum_wavelength_for_correction
        wavelength_max = self.maximum_wavelength_for_correction
        a = 0
        b = 0
        i = 0  # set the starting counter
        z = 0  # set the ending counter

        """This loop set the steps between wavelength_min & max."""
        for element in table_wavelength.flat:
            if element <= float(wavelength_max):
                z += 1
            if element < float(wavelength_min):
                i += 1

        """Create a table of the sample containing only the intensity between wavelength min & max."""
        emission_limited = table_fluo[i:z]

        """Check the length of the array emission_limited and corr_fact, then correct the USB spectro response,
        then correct the area of the emission."""
        for _ in emission_limited:
            a += 1
        for _ in corr_fact:
            b += 1
        if a == b:
            emission_limited_ccd_corr = emission_limited * corr_fact
            area_fluo = np.trapz(emission_limited_ccd_corr)
            area_fluo_corr = area_fluo * corr_dye_area

            """Return area_fluo_corr: float
            This float will be used for sigma 2 calculation"""

            logging.info("Tpa_calc:_Emission_Correction => ok !")

            return area_fluo_corr

    def _Quad_Log(self, list_of_fully_corrected_area, list_of_square_dark_corr_power):
        """
        This function determine the quadraticity of the F/P²
        :param list_of_fully_corrected_area: lst
        :param list_of_square_dark_corr_power: lst
        :return: int, float, float, float, array, array, array
        """

        logging.info("Tpa_calc:Starting _Quad_Log.")
        logging.info(f"Tpa_calc:list_of_fully_corrected_area: {list_of_fully_corrected_area}.")
        logging.info(f"Tpa_calc:list_of_square_dark_corr_power: {list_of_square_dark_corr_power}.")

        try:
            reg_lin= linregress(list_of_square_dark_corr_power, list_of_fully_corrected_area)
            area= round(reg_lin[0],1)
            correlation_coeff=round(reg_lin[2]*reg_lin[2],5)
            fluo_log = []
            power_log = []
            for fluo in list_of_fully_corrected_area:

                logging.info(f"Tpa_calc:log fluo: {np.log(fluo)}.")

                fluo_log.append(np.log(fluo))
            for power in list_of_square_dark_corr_power:

                logging.info(f"Tpa_calc:log power= {np.log(power)}.")
                logging.info(f"Tpa_calc:log (F)/log(P²)= {np.log(fluo) / np.log(power)}.")

                power_log.append(np.log(power))
            fluo_log= np.array(fluo_log)
            power_log = np.array(power_log)
            slope, intercept, pearson, pvalue, std_dev = linregress(power_log, fluo_log)
            coeff = pearson*pearson  # Here I recover r² of the fitting
            fit = intercept + slope*power_log
        except:

            logging.error("Tpa_calc:Error.")

            return 0, 0, 0, 0, [], [], 0

        logging.info("Tpa_calc:_Quad_Log => ok !")

        return area, correlation_coeff, slope, coeff, fluo_log, power_log, fit

    def _Figure_Duo(self, x1, y1, x2, y2, name):
        """
        This function draws figures containing 2 plots (x1,y1) and (x2,y2) and save them
        as name.png to root/Figures folder.
        :param x1: array
        :param y1: array
        :param x2: array
        :param y2: array
        :param name: str
        :return: bool
        """

        logging.info("Starting _Figure_Duo.")
        logging.info(f"Tpa_calc:type(x1): {type(x1)}.")
        logging.info(f"Tpa_calc:type(y1): {type(y1)}.")
        logging.info(f"Tpa_calc:type(x2): {type(x2)}.")
        logging.info(f"Tpa_calc:type(y2): {type(y2)}.")
        logging.info(f"Tpa_calc:name : {name}.")

        fig, ax = plt.subplots()
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter(
            '%.1e'))  # here %.1e is used to put scientific notation  with 1 significative value
        ax.xaxis.set_major_formatter(
            ticker.FormatStrFormatter('%.0f'))  # here %.0f is used to put decimal with 0 significative value
        ax.xaxis.set_major_locator(ticker.MaxNLocator(10))
        ax.xaxis.set_minor_locator(ticker.MaxNLocator(100))
        ax.yaxis.set_major_locator(ticker.MaxNLocator(6))
        ax.plot(x1, y1, color='black')
        ax.plot(x2, y2, color='red')
        plt.xlabel('Wavelength / nm')
        plt.ylabel(r'Fluorescence intensity / cps')
        with open(ROOT,"r") as r:
            root=json.load(r)
        root=Path(f"{root}")
        root = root / "Datas"/ "Figures"
        root.mkdir(exist_ok=True)
        fig.savefig(root / f"{name}.png")

        logging.info(f"Tpa_calc:Saving figure in {root} folder.")
        logging.info("Tpa_calc:Figure created and saved => ok !")

        return True
