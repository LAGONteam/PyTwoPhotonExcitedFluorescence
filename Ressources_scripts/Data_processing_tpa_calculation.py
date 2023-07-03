"""
Created on 5/1/2023
@author: Jonathan DANIEL, jonathan.daniel@u-bordeaux.fr
"""
import copy
from pathlib import Path
import os.path
import logging
import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
from scipy.stats import linregress
import csv
import time

from Ressources_scripts import References_data
from Ressources_scripts.Data_processing import save

import logging

logging.basicConfig(level=logging.INFO,
                    filename="PyTPEF.log",
                    filemode="a",
                    format='%(asctime)s - %(levelname)s - %(message)s')

CUR_DIR = os.path.dirname(__file__)
DATA = Path(f"{CUR_DIR}").parent
root= DATA / "Datas"
root.mkdir(exist_ok=True)

class Read_Data_File():

    def __init__(self):
        """Dictionaries"""
        self.refractive_index_dic = {"Water": 1.333, "Toluene": 1.4961, "THF": 1.405, "Methanol": 1.3288, "DMSO": 1.477, "DCM": 1.424, "Ethanol":1.361, "DMF":1.429, "Acetonitrile": 1.344, "Acetone": 1.359, "1,2-Dichloroethane": 1.442, "CHCl3": 1.444, "Acetic acid": 1.327, "1,4-Dioxane": 1.422, "Cyclohexane": 1.427, "Heptane":1.394}
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

        """Variables"""
        self.root=root

    def _Calcul_TPA(self, sample_info, processed_data, sample_name, reference_name, wavelength, simulation, root):
        """
        This function calculates the two-photon cross section sigma 2 of sample at a given wavelength.
        :param sample_info: dict
        :param processed_data: dict
        :param sample_name: str
        :param reference_name: str
        :param wavelength: int
        :param simulation: bool
        :param root: path
        :return: float, float
        """

        logging.info("RDF:Starting _Calcul_TPA.")
        logging.info(f"RDF:sample_info: {sample_info}.")
        logging.info(f"RDF:processed_data: {processed_data}.")
        logging.info(f"RDF:sample_name: {sample_name}.")
        logging.info(f"RDF:reference_name: {reference_name}.")
        logging.info(f"RDF:wavelength: {wavelength}.")
        logging.info(f"RDF:simulation: {simulation}.")
        logging.info(f"RDF:root: {root}.")
        if sample_name == reference_name:
            ref_data= processed_data[sample_name]
            logging.info("RDF:Reference and sample are the same !")
            logging.info(f"RDF:sample_name: {sample_name}, reference_name: {reference_name}.")
            logging.info(f"RDF:type(sample_name): {type(sample_name)}, type(reference_name): {type(reference_name)}.")
        else:
            try:
                ref_data= save().read_process(sample= reference_name, simulation=simulation, root_=root)
                logging.info(f"RDF:Ref_data: {ref_data}.")
            except FileNotFoundError:
                logging.error("RDF:File not Found. Please measure the reference before the sample !")
                print("ERROR_1")
                return -1, -1

        try:
            logging.info(f"RDF:sample_name: {sample_name}.")
            logging.info(f"RDF:processed_data.keys(): {processed_data.keys()}.")
            logging.info(f"RDF:processed_data[sample_name]: {processed_data[sample_name]}.")
            logging.info(f"RDF:processed_data[reference_name]: {processed_data[reference_name]}.")
        except KeyError:
            print("ERROR_2")
            logging.warning("RDF:Key error")

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
            area = (processed_data[f"{sample_name}"][f"{wavelength}"]['full_corrected_area'])/(ref_data[f"{wavelength}"]['full_corrected_area'])
            logging.info(f"RDF:area: {area}.")
        except KeyError:
            logging.warning("RDF:Key Error. Please measure the reference before the sample !")
            print("ERROR_3")
            return -1, -1
        concentration = float(ref_concentration) / float(sample_concentration)
        refractive_index = reference_refractive_index / sample_refractive_index
        calcul_1 = area * concentration * refractive_index
        logging.info(f"RDF:reference_sigma2.keys(): {reference_sigma2.keys()}.")
        time.sleep(0.1)
        calcul_tpa = calcul_1 * (reference_sigma2[f"{wavelength}"])
        logging.info(f"RDF:calculation 2PA: {calcul_tpa}, type (calcul_tpa): {type(calcul_tpa)}.")
        logging.info(f"RDF:Phi: {phi}, type(phi): {type(phi)}.")
        calcul_s2 = calcul_tpa / phi
        sigma2_phi = round(calcul_tpa, 1)
        sigma2 = round(calcul_s2, 1)
        logging.info("RDF:_Calcul_TPA => ok !")

        area_ref= ref_data[f"{wavelength}"]['full_corrected_area']
        area_sample=processed_data[f"{sample_name}"][f"{wavelength}"]['full_corrected_area']
        s2_ref_=reference_sigma2[f"{wavelength}"]
        print("/"*1000)
        print(
            "*"*1000,"\n",
            "wavelength", wavelength,"\n",
            "sample", sample_name,"\n",
            "ref_area", area_ref,"\n",
            "area_sample",area_sample,"\n",
            "area", area,"\n",
            "ref_conc", ref_concentration,"\n",
            "ref_sample", sample_concentration,"\n",
            "concentration", concentration,"\n",
            "ref_index", reference_refractive_index,"\n",
            "sample_index", sample_refractive_index,"\n",
            "refractive index", refractive_index,"\n",
            "calcul_1", calcul_1,"\n",
            "s2_ref_", s2_ref_,"\n",
            "calcul_tpa", calcul_tpa,"\n",
            "calcul_s2", calcul_s2,"\n",
            "*" * 1000
            )

        return sigma2_phi, sigma2

    def reference_dye_info(self, ref_dye):
        """
        unused function
        """
        """
        This function returns informations about the reference dye.
        :param ref_dye: str
        :return: dict, str, float
        """
        logging.info("RDF:Starting reference_dye_info")
        logging.info(f"RDF:ref_dye: {ref_dye}.")
        if ref_dye == "fluorescein":
            solvent = "Water"
            phi = 0.9
            logging.info("reference_dye_info => ok !")
            return self.spectrum_fluo_fluorescein_water, solvent, phi
        elif ref_dye == "NR":
            solvent = "DMSO"
            phi = 0.79
            logging.info("RDF:reference_dye_info => ok !")
            return self.spectrum_fluo_NR_DMSO, solvent, phi

    def _Read_Reference_Fluo(self, dye_for_correction, dye_as_reference):
        """
        unused function
        """
        """
        This function reads and calculate the corrected area factor (depending on the CCD dye)
        :param dye_for_correction: str
        :param dye_as_reference: str
        :return: float
        """
        logging.info("RDF:Starting _Read_Reference_Fluo")
        logging.info(f"RDF:dye_for_correction: {dye_for_correction}.")
        logging.info(f"RDF:dye_as_reference: {dye_as_reference}.")

        ref = dye_as_reference
        if dye_for_correction == "JOD":
            wavelength_min = 480.0
            wavelength_max = 650.0
        elif dye_for_correction == "DMANs":
            wavelength_min = 510.0
            wavelength_max = 750.0

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
        logging.info(f"RDF:area_full_dye: {area_full_dye}.")
        logging.info(f"RDF:area_limited_dye: {area_limited_dye}.")
        corr_dye_area = area_full_dye / area_limited_dye
        logging.info(f"RDF:corr_dye_area: {corr_dye_area}.")

        self._Figure_Duo(self.dye_wavelength, self.dye_fluo, self.dye_limited_wavelength, self.dye_limited_fluo, f"{ref} limited")
        self.dye_fluo.clear()
        self.dye_limited_fluo.clear()
        self.dye_limited_wavelength.clear()
        self.dye_wavelength.clear()
        logging.info("RDF:_Read_Reference_Fluo => ok !")
        return corr_dye_area

    def _Read_Dyes_Emission_Spectrum(self, root, dye_for_correction, sample_name):
        """
        unused function
        """
        """
        This function reads the fluorescence spectrum of the sample and calculate the corrected area factor
        (depending on the CCD dye)
        :param root: path
        :param dye_for_correction: str
        :param sample_name: str
        :return: float
        """

        logging.info("RDF:Starting _Read_Dyes_Emission_Spectrum.")
        logging.info(f"RDF:root: {root}")
        logging.info(f"RDF:dye_for_correction: {dye_for_correction}")
        logging.info(f"RDF:sample_name: {sample_name}")
        if dye_for_correction == "JOD":
            wavelength_min = 480.0
            wavelength_max = 650.0
        elif dye_for_correction == "DMANs":
            wavelength_min = 510.0
            wavelength_max = 750.0

        logging.info(f"RDF:dye_for_correction:{dye_for_correction}.")
        logging.info(f"RDF:wavelength_min: {wavelength_min}.")
        logging.info(f"RDF:wavelength_max: {wavelength_max}.")

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
        logging.info(f"RDF:area_full_dye: {area_full_dye}")
        logging.info(f"RDF:area_limited_dye: {area_limited_dye}")


        """Step 3. Calculate the area correction factor to takes into account the part of the fluorescence intensity
        which is not take into account by the ccd correction (out of range)."""
        corr_dye_area = area_full_dye / area_limited_dye
        logging.info(f"RDF:corr_dye_area: {corr_dye_area}")

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
        logging.info("RDF:_Read_Dyes_Emission_Spectrum => ok !")
        return corr_dye_area

    def _Spectro_Correction_Factor(self, dye_for_correction, x, y_corr):
        """
        unused function
        """
        """
        This function corrects the spectral response of the USB spectro using a fluorescent dye. The wavelength range
        (min & max) depends on the emission spectrum of the fluorescent dye => no correction where the dye is not
        fluorescent.
        :param dye_for_correction: str
        :param x: array
        :param y_corr: array
        :return: lst, float
        """

        logging.info("RDF:Starting _Spectro_Correction_Factor.")
        logging.info(f"RDF:dye_for_correction: {dye_for_correction}.")
        logging.info(f"RDF:x: {x}.")
        logging.info(f"RDF:y_corr: {y_corr}.")
        if dye_for_correction == "JOD":
            ref_dye = self.spectrum_JOD_water
            wavelength_min = 480
            wavelength_max = 650
        elif dye_for_correction == "DMANs":
            ref_dye = self.spectrum_DMANs_toluen
            wavelength_min = 510
            wavelength_max = 750
        wavelength_reference_tpa = []
        fluo_reference_tpa = []

        """Loop to check that the arrays will have the same size with same wavelengths, require that the emission 
        spectrum of the ccd correction dye is record with a pitch of 0.1 nm."""
        for _ in x.flat:
            _ = round(float(_), 1)
            if _ > wavelength_min and _ < wavelength_max:
                wavelength_reference_tpa.append(f"{_}")
                fluo_reference_tpa.append(ref_dye[f"{_}"])
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
        logging.info(f"RDF:wavelength_measured_sliced.size: {wavelength_measured_sliced.size}, wavelength_measured_sliced: {wavelength_measured_sliced}.")
        logging.info(f"RDF:len(wavelength_reference_tpa): {len(wavelength_reference_tpa)}, wavelength_reference_tpa: {wavelength_reference_tpa}.")

        if wavelength_measured_sliced.size == len(wavelength_reference_tpa):
            logging.info(f"RDF:fluroescence intensity sliced: {intensity_measured_sliced}.")
            corr_fact = intensity_ref / intensity_measured_sliced
            logging.info(f"RDF:corrector factor minimum is: {corr_fact.min()}.")
            corr_fact /= (corr_fact.min())

        else:
            logging.error("RDF:Error in wavelength_measured_sliced.size")
            corr_fact = 0
        """Return corr_fact: array 
        This array correspond to a correction of spectro response between wavelength_min & max only !"""
        logging.info("RDF:_Spectro_Correction_Factor => ok !")
        return wavelength_measured_sliced, corr_fact

    def _Emission_Correction(self, corr_dye_area, table_fluo, table_wavelength, dye_for_correction, corr_fact):
        """
        unused function
        """
        """
        Correct the emission spectrum of the sample with corr_fact (USB spectro response), then correct the emission
        area, and return a fully corrected intensity area for sigma 2 calculation)

        :param corr_dye_area: float
        :param table_fluo: array
        :param table_wavelength: array
        :param dye_for_correction: str
        :param corr_fact: array
        :return: float
        """

        logging.info("RDF: Starting _Emission_Correction.")
        logging.info(f"RDF:corr_dye_area: {corr_dye_area}.")
        logging.info(f"RDF:table_fluo: {table_fluo}.")
        logging.info(f"RDF:table_wavelength: {table_wavelength}.")
        logging.info(f"RDF:dye_for_correction: {dye_for_correction}.")
        logging.info(f"RDF:corr_fact: {corr_fact}.")

        if dye_for_correction == "JOD":
            wavelength_min = 480
            wavelength_max = 650
        elif dye_for_correction == "DMANs":
            wavelength_min = 510
            wavelength_max = 750

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
            logging.info("RDF:_Emission_Correction => ok !")
            return area_fluo_corr

    def _Quad_Log(self, list_of_fully_corrected_area, list_of_square_dark_corr_power):
        """
        unused function
        """

        logging.info("RDF:Starting _Quad_Log.")
        logging.info(f"RDF:list_of_fully_corrected_area: {list_of_fully_corrected_area}.")
        logging.info(f"RDF:list_of_square_dark_corr_power: {list_of_square_dark_corr_power}.")
        reg_lin= linregress(list_of_square_dark_corr_power, list_of_fully_corrected_area)
        area= round(reg_lin[0],1)
        correlation_coeff=round(reg_lin[2]*reg_lin[2],5)

        temp_list=[]

        for n in range(len(list_of_fully_corrected_area)):
            print(n)
            temp_area =list_of_fully_corrected_area[n]/list_of_square_dark_corr_power[n]
            temp_list.append(temp_area)
        print(temp_list)

        area_bis = np.mean(temp_list)
        print(area_bis)
        print(area)

        fluo_log = []
        power_log = []
        for fluo in list_of_fully_corrected_area:
            logging.info(f"RDF:log fluo: {np.log(fluo)}.")
            fluo_log.append(np.log(fluo))
        for power in list_of_square_dark_corr_power:
            logging.info(f"RDF:log power: {np.log(power)}.")
            logging.info(f"RDF:log (F)/log(P²): {np.log(fluo)/np.log(power)}.")
            power_log.append(np.log(power))
        fluo_log= np.array(fluo_log)
        power_log = np.array(power_log)
        slope, intercept, pearson, pvalue, std_dev = linregress(power_log, fluo_log)
        coeff = pearson*pearson  # Here I recover r² of the fitting
        fit = intercept + slope*power_log
        logging.info("RDF:_Quad_Log => ok !")

        return area_bis, correlation_coeff, slope, coeff, fluo_log, power_log, fit

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
        logging.info(f"RDF:x1: {x1}.")
        logging.info(f"RDF:y1: {y1}.")
        logging.info(f"RDF:x2: {x2}.")
        logging.info(f"RDF:y2: {y2}.")
        logging.info(f"RDF:name : {name}.")
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
        root = self.root / "Figures"
        root.mkdir(exist_ok=True)
        logging.info(f"RDF:Saving figure in {root} folder.")
        fig.savefig(root / f"{name}.png")
        logging.info("RDF:_Figure_Duo => ok !")
        return True