"""
Created on 5/1/2023
@author: Jonathan DANIEL, jonathan.daniel@u-bordeaux.fr
"""
import copy
import os
import json
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import ticker
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import linregress
import logging

logging.basicConfig(level=logging.INFO,
                    filename="PyTPEF.log",
                    filemode="a",
                    format='%(asctime)s - %(levelname)s - %(message)s')

CUR_DIR = os.path.dirname(__file__)
DATA = Path(f"{CUR_DIR}").parent
DATA= DATA / "Datas"
DATA.mkdir(exist_ok=True)
root = DATA / "Metadata" / "sample_data.json"

def get_samples(root):
    """
    Unused function
    """

    logging.info("Data_proc:Starting get_samples.")
    logging.info(f"Data_proc:root: {root}.")

    samples = []
    with open(root, "r") as f:
        sample_names = json.load(f)
        for _ in sample_names:
            samples.append(sample_data(_))

        logging.info("Data_proc:get_samples => ok !")

        return samples

def get_all_samples():
    """
    Unused function
    """

    logging.info("Data_proc:Starting get_all_samples.")

    samples = []
    with open(root, "r") as f:
        sample_names = json.load(f)
        for _ in sample_names:
            samples.append(sample_data(_))

        logging.info("Data_proc:get_all_samples => ok !")

        return samples

def get_wavelength(sample):
    """
    Unused function
    """

    logging.info("Data_proc:Starting get_wavelength.")
    logging.info(f"Data_proc:sample: {sample}.")

    wavelength = []
    root_in = DATA / f"{sample}" / "wavelength.json"
    with open(root_in, "r") as f:
        wavelengths = json.load(f)
        for _ in wavelengths:
            wavelength.append(sample_data(str(_)))

        logging.info("Data_proc:get_wavelength => ok !")

        return wavelength

class Temp_Files():

    def __init__(self):
        pass

    def save_files(self, data, name):
        """
        Unused function
        """

        logging.info("Data_proc:Starting save_files.")
        logging.info(f"Data_proc:data: {data}.")
        logging.info(f"Data_proc:name: {name}.")

        root_in = DATA / "temp"
        root_in.mkdir(exist_ok=True)
        root_in = root_in / f"{name}.json"
        root_in.touch()
        if name =="root_dict":
            data= f"{data}"
        with open(root_in, "w") as w:
            json.dump(data, w, cls= NumpyEncoder)

        logging.info("Data_proc:save_files => ok !")

    def get_files(self, name):
        """
        Unused function
        """

        logging.info("Data_proc:Starting get_files.")
        logging.info(f"Data_proc:name: {name}.")

        root_in = DATA / "temp" / f"{name}.json"
        with open(root_in, "r") as r:
            data = json.load(r)

        logging.info("Data_proc:get_files => ok !")

        return data

class sample_data:

    def __init__(self, sample_name):
        self.root = root
        self.sample_name = sample_name
        with open(self.root, 'w') as f:
            json.dump([], f)

    def __str__(self):
        return self.sample_name

    def _get_sample_data(self):
        """
        This function returns the content of the sample_data.json file
        :return: dict
        """
        with open(self.root, "r") as f:
            return json.load(f)

    def _write_sample_data(self, sample_name):
        """
        This function writes the data in the sample_data.json file
        :param sample_name: str
        :return: None
        """
        with open(self.root, 'w') as f:
            json.dump(sample_name, f)

    def add_to_sample_name(self):
        """
        Unused function
        """
        data = self._get_sample_data()
        data.append(self.sample_name)
        self._write_sample_data(data)
        return True

    def remove_sample(self):
        """
        This function removes the name of the sample from the sample_data.json file
        :return: None
        """
        data = self._get_sample_data()

        if self.sample_name in data:
            data.remove(self.sample_name)
            self._write_sample_data(data)

class NumpyEncoder(json.JSONEncoder):
    """ Special json encoder for numpy types """
    def default(self, obj):
        if isinstance(obj, (np.int_, np.intc, np.intp, np.int8,
                            np.int16, np.int32, np.int64, np.uint8,
                            np.uint16, np.uint32, np.uint64)):
            return int(obj)
        elif isinstance(obj, (np.float_, np.float16, np.float32,
                              np.float64)):
            return float(obj)
        elif isinstance(obj, (np.ndarray,)):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)

class save:

    def __init__(self):
        pass

    def read_column_from_excel_files(self, sample, root):
        """
        Unused function
        """

        logging.info("Data_proc:Starting read_column_from_excel_files.")
        logging.info(f"Data_proc:sample: {sample}.")
        logging.info(f"Data_proc:root: {root}.")

        file_name=root / f"{sample} informations.xlsx"
        df= pd.read_excel(io=file_name)
        intensity= df['Reference intensity / cps']
        wavelength= df['Reference wavelength / nm']
        intensity= np.array(intensity)
        wavelength= np.array(wavelength)

        logging.info("Data_proc:read_column_from_excel_files => ok !")

        return wavelength, intensity

    def _get_data_for_new_mesure(self, sample, wavelength, simulation, root):
        """
        This function gives the data of the sample at a given wavelength.
        :param sample: str
        :param wavelength: int
        :param simulation: bool
        :param root: path
        :return: dict
        """

        logging.info("Data_proc:Starting _get_data_for_new_mesure.")
        logging.info(f"Data_proc:sample: {sample}.")
        logging.info(f"Data_proc:wavelength: {wavelength}.")
        logging.info(f"Data_proc:simulation: {simulation}.")
        logging.info(f"Data_proc:root: {root}.")

        if simulation == False:
            root = root / f"{sample}" / f"experimental_data {sample}_{wavelength}.json"
        elif simulation == True:
            root = root / f"{sample}" / f"experimental_data {sample}_{wavelength}_temp.json"
        with open(root, "r") as r:
            data = json.load(r)

        logging.info("Data_proc:_get_data_for_new_mesure => ok !")

        return data

    def _get_file_for_new_mesure(self, sample, wavelength, root):
        """
        This function give the path for a new measure of the sample at a given wavelength.
        :param sample: str
        :param wavelength: int
        :param root: path
        :return: path
        """

        logging.info("Data_proc:Starting _get_file_for_new_mesure.")
        logging.info(f"Data_proc:sample: {sample}.")
        logging.info(f"Data_proc:wavelength: {wavelength}.")
        logging.info(f"Data_proc:root: {root}.")

        root_in = root / f"{sample}" / f"experimental_data {sample}_{wavelength}.json"

        logging.info("Data_proc:_get_file_for_new_mesure => ok !")

        return root_in

    def fit_curve(self, x, y):
        """
        This function fits the curve log(y) =f(log(x²)) using a 2nd order polynome.
        The slope corresponds to the order of the absorption process
        :param x: array
        :param y: array
        :return: array, array, str, float, array, array, array, array
        """
        try :

            logging.info("Data_proc:Starting fit_curve.")
            logging.info(f"Data_proc:type(x): {type(x)}.")
            logging.info(f"Data_proc:type(y): {type(y)}.")

            line, stats= np.polynomial.polynomial.polyfit(x=x, y=y, deg=2, full=True)
            #line = np.poly1d(np.polyfit(x, y, 1))

            logging.info(f"Data_proc:line: {line}.")
            logging.info(f"Data_proc:stats: {stats}.")

            x_expand= np.linspace(min(x), max(x), 100)
            y_expand = line[0] + line[1]*x_expand + line[2]*x_expand*x_expand

            x_log= np.log(x)
            y_log= np.log(y)

            logar, stats_log=np.polynomial.polynomial.polyfit(x= x_log, y=y_log, deg=1, full=True)

            logging.info(f"Data_proc:logar: {logar}.")
            logging.info(f"Data_proc:stats_log: {stats_log}.")

            order_of_absorption= round(np.abs(logar[1]),2)

            x_log2= np.linspace(np.log(min(x)), np.log(max(x)), 100)
            y_log2= logar[0] + logar[1]*x_log2

            # logar[1] absorption order
            a=round(line[0],0)
            b=round(line[1],1)
            c=round(line[2],1)
            equation = f"{a:.1E} + {b:.1E}*x + {c:.1E}*x²" #fstring :1.E allow scientific notation with 1 digit

            logging.info("Data_proc:fit_curve => ok !")

        except :

            logging.warning("Data_proc:fit_curve => Problem !")

            return x, y, "Error", 0, x, y, x, y


        return x_expand, y_expand, equation, order_of_absorption, x_log, y_log, x_log2, y_log2

    def read_samples_informations(self, root):
        """
        This function reads the data stored in the samples and measure informations.json file
        :param root: path
        :return: dict
        """

        logging.info("Data_proc:Starting read_samples_informations.")
        logging.info(f"Data_proc:root: {root}.")

        root =root / "Figures" / "Experimental data" / "samples and measure informations.json"
        root.touch()
        with open(root, "r") as r:
            data= json.load(r)

        logging.info("Data_proc:read_samples_informations => ok !")

        return data

    def change_samples_informations(self, root, new_data):
        """
        Unused function
        """

        logging.info("Data_proc:Starting change_samples_informations.")
        logging.info(f"Data_proc:root: {root}.")
        logging.info(f"Data_proc:type(new_data): {type(new_data)}.")

        root = root / "Figures" / "Experimental data" / "samples and measure informations.json"
        with open(root, "w") as w:
            json.dump(new_data, w)

        logging.info("Data_proc:change_samples_informations => ok !")

        return True

    def update_samples_informations(self, original_data, new_data, type_of_data, sample, root, new_data_name=None):
        """
        This function is not ready yet, it will be used to allow user to change some parameters data after measurement.
        """

        logging.info("Data_proc:Starting update_samples_informations.")
        logging.info(f"Data_proc:original_data: {original_data}.")
        logging.info(f"Data_proc:new_data: {new_data}.")
        logging.info(f"Data_proc:type_of_data: {type_of_data}.")
        logging.info(f"Data_proc:sample: {sample}.")

        data = self.read_samples_informations(root=root)

        logging.info(f"Data_proc:data before update: {data}.")

        for key in data.keys():
            if key == str(sample):
                for key_ in data[f"{sample}"].keys():
                    if key_ == type_of_data:
                        if new_data_name != None:
                            data[f"{sample}"][f"{key_}"]={new_data_name:new_data}
                        else:
                            data[f"{sample}"][f"{key_}"]=new_data

        logging.info(f"Data_proc:data after update: {data}.")

        self.change_samples_informations(root=root, new_data=data)

        logging.info("Data_proc:update_samples_informations => ok !")


    def copy_data_for_simulation(self, sample, wavelength, root):
        """
        This functions creates a temporary file of experimental_data sample wavelength.json to store
        simulation data
        :param sample: str
        :param wavelength: int
        :param root: path
        :return: None
        """

        logging.info("Data_proc:Starting copy_data_for_simulation.")
        logging.info(f"Data_proc:sample: {sample}.")
        logging.info(f"Data_proc:wavelength: {wavelength}.")
        logging.info(f"Data_proc:root: {root}.")

        root_in = root / f"{sample}" / f"experimental_data {sample}_{wavelength}.json"
        with open(root_in, "r") as r:
            data = json.load(r)
        root_out = root / f"{sample}" / f"experimental_data {sample}_{wavelength}_temp.json"
        root_out.touch()
        with open(root_out, "w") as w:
            json.dump(data, w)

        logging.info("Data_proc:copy_data_for_simulation => ok !")

    def copy_processed_data_for_simulation(self, sample, root):
        """
        This functions creates a temporary file of processed_data sample wavelength.json to store
        simulation data
        :param sample: str
        :param root: path
        :return: None
        """

        logging.info("Data_proc:Starting copy_processed_data_for_simulation.")
        logging.info(f"Data_proc:sample: {sample}.")
        logging.info(f"Data_proc:root: {root}.")

        root_in = root / f"{sample}" / f"processed_data {sample}.json"
        with open(root_in, "r") as r:
            data = json.load(r)
        root_out = root / f"{sample}" / f"processed_data {sample}_temp.json"
        root_out.touch()
        with open(root_out, "w") as w:
            json.dump(data, w)

        logging.info("Data_proc:copy_processed_data_for_simulation => ok !")

    def extract_data_from_experimental(self, sample, wavelength, root):
        """
        This function extracts data from the experimental_data sample wavelength.json file.
        :param sample: str
        :param wavelength: int
        :param root: path
        :return: dict
        """

        logging.info("Data_proc:Starting extract_data_from_experimental.")
        logging.info(f"Data_proc:sample: {sample}.")
        logging.info(f"Data_proc:wavelength: {wavelength}.")
        logging.info(f"Data_proc:root: {root}.")

        root = root / f"{sample}" / f"experimental_data {sample}_{wavelength}.json"
        with open(root, "r") as r:
            data = json.load(r)

        logging.info("Data_proc:extract_data_from_experimental => ok !")

        return data

    def extract_data_from_experimental_bis(self,sample, wavelength, root):
        """
        This function is the same as below, check why !
        :param sample:
        :param wavelength:
        :param root:
        :return:
        """

        logging.info("Data_proc:Starting extract_data_from_experimental_bis.")
        logging.info(f"Data_proc:sample: {sample}.")
        logging.info(f"Data_proc:wavelength: {wavelength}.")
        logging.info(f"Data_proc:root: {root}.")

        root= root / f"{sample}" / f"experimental_data {sample}_{wavelength}.json"
        with open(root, "r") as r:
            data = json.load(r)

        logging.info("Data_proc:extract_data_from_experimental_bis => ok !")

        return data

    def _get_intensity_vs_wavelength_data(self, data):
        """
        This functions extracts wavelength and intensity-dark arrays from data
        :param data: dict
        :return: array, array
        """

        logging.info("Data_proc:Starting _get_intensity_vs_wavelength_data.")
        logging.info(f"Data_proc:type(data): {type(data)}.")

        x = data['wavelength']
        y = data['intensity - dark']

        logging.info("Data_proc:_get_intensity_vs_wavelength_data => ok !")

        return x, y

    def _get_quick_f_and_pp(self, sample, wavelength, root, number_of_measure):
        """
        This functions extracts intensity and square power arrays from experimental data
        of sample at a specific wavelength.
        :param sample: str
        :param wavelength: int
        :param root: path
        :param number_of_measure: int
        :return: lst, lst
        """

        logging.info("Data_proc:Starting _get_quick_f_and_pp.")
        logging.info(f"Data_proc:sample: {sample}.")
        logging.info(f"Data_proc:wavelength: {wavelength}.")
        logging.info(f"Data_proc:root: {root}.")
        logging.info(f"Data_proc:number_of_measure: {number_of_measure}.")

        data_bis= self.extract_data_from_experimental(sample=sample, wavelength=wavelength, root=root)
        list_of_fluo=[]
        list_of_pp=[]
        for n in range(number_of_measure):
            try:
                temp_var = data_bis[f"{wavelength}"].keys()

                logging.info(f"Data_proc:type(data_bis): {type(data_bis)}.")
                logging.info(f"Data_proc:data_bis[wavelength].keys(): {temp_var}.")

                list_of_fluo.append(data_bis[f"{wavelength}"][f"{n}"]["fully corrected area"])
                list_of_pp.append(data_bis[f"{wavelength}"][f"{n}"]["raw power - dark"])

            except KeyError:

                logging.warning("Data_proc:Invalid value of power pitch.")

        logging.info("Data_proc:_get_quick_f_and_pp => ok !")

        return list_of_fluo, list_of_pp

    def extract_process_data(self, sample, wavelength, simulation, root):
        """
        This functions extracts datas from processed data of sample
        :param sample: str
        :param wavelength: int
        :param simulation: bool
        :param root: path
        :return: lst, lst, lst, lst, lst, lst
        """

        logging.info("Data_proc:Starting extract_process_data.")
        logging.info(f"Data_proc:sample: {sample}.")
        logging.info(f"Data_proc:wavelength: {wavelength}.")
        logging.info(f"Data_proc:simulation: {simulation}.")
        logging.info(f"Data_proc:root: {root}.")

        list_of_s2f = []
        list_of_wavelength = []
        if simulation == False:
            root_in = root / f"{sample}" / f"processed_data {sample}.json"
        elif simulation == True:
            root_in = root / f"{sample}" / f"processed_data {sample}_temp.json"
        with open(root_in, "r") as r:
            data = json.load(r)

            logging.info(f"Data_proc:extract_process_data: {data.keys()}.")

        for key_wavelength in data:
            temp_var = data[f"{key_wavelength}"].keys()

            logging.info(f"Data_proc:extract_process_data_bis: {temp_var}.")

            list_of_s2f.append(data[f"{key_wavelength}"]["S2F"])
            list_of_wavelength.append(key_wavelength)
        if wavelength == None:
            list_of_log_f = False
            list_of_log_pp = False
            list_of_fluo = False
            list_of_pp = False
        else:
            list_of_log_f= data[f"{wavelength}"]["quadraticity"]["log(F)"]
            list_of_log_pp= data[f"{wavelength}"]["quadraticity"]["log(P\u00b2)"]
            list_of_fluo=data[f"{wavelength}"]["quadraticity"]['F']
            list_of_pp= data[f"{wavelength}"]["quadraticity"]['PP']

        logging.info("Data_proc:extract_process_data => ok !")

        return list_of_s2f, list_of_wavelength, list_of_log_f, list_of_log_pp, list_of_fluo, list_of_pp

    def remove_data(self, sample, wavelength, measure, simulation, root):
        """
        This functions removes data from experimental_data sample wavelength (temp).json file
        :param sample: str
        :param wavelength: int
        :param measure: int
        :param simulation: bool
        :param root: path
        :return: None
        """

        logging.info("Data_proc:Starting remove_data.")
        logging.info(f"Data_proc:sample: {sample}.")
        logging.info(f"Data_proc:wavelength: {wavelength}.")
        logging.info(f"Data_proc:measure: {measure}.")
        logging.info(f"Data_proc:simulation: {simulation}.")
        logging.info(f"Data_proc:root: {root}.")

        if simulation == False:
            root_in = root / f"{sample}" / f"experimental_data {sample}_{wavelength}.json"
        elif simulation == True:
            root_in = root / f"{sample}" / f"experimental_data {sample}_{wavelength}_temp.json"
        with open(root_in, "r") as r:
            data = json.load(r)
        data[f"{wavelength}"].pop(f"{measure}")

        logging.info("Data_proc:Measure delete succesfuly !")

        with open(root_in, "w") as w:
            json.dump(data, w)

        logging.info("Data_proc:remove_data => ok !")

    def _get_location_file(self, root):
        """
        Unused function
        """

        logging.info(f"Data_proc:Starting _get_location_file: root: {root}.")

        return root

    def write_only(self, dictionary_process, sample, root):
        """
        Unused function
        """

        logging.info("Data_proc:Starting write_only.")
        logging.info(f"Data_proc:type(dictionary_process): {type(dictionary_process)}.")
        logging.info(f"Data_proc:sample: {sample}.")
        logging.info(f"Data_proc:root: {root}.")

        root_in = root / sample
        root_in.mkdir(exist_ok=True)
        root = root_in / f"processed_data {sample}.json"
        root.touch()
        with open(root, "w") as w:
            json.dump(dictionary_process, w, cls=NumpyEncoder)

        logging.info("Data_proc:write_only => ok !")

    def extract_datas(self, dictionnary):
        """
        Unused function
        """

        logging.info("Data_proc:Starting extract_datas.")
        logging.info(f"Data_proc:type(dictionnary): {type(dictionnary)}.")

        x = []
        y = []
        for key in dictionnary:
            value = dictionnary[key]["quadraticity"]["slope"]
            x.append(key)
            y.append(value)
        x = np.array(x)
        y = np.array(y)

        logging.info("Data_proc:extract_datas => ok !")

        return x, y

    def save_exp(self, dictionary, sample, wavelength, root):
        """
        Unused function
        """

        logging.info("Data_proc:Starting save_exp.")
        logging.info(f"Data_proc:type(dictionary): {type(dictionary)}.")
        logging.info(f"Data_proc:sample: {sample}.")
        logging.info(f"Data_proc:wavelength: {wavelength}.")
        logging.info(f"Data_proc:root: {root}.")

        root_in = root / sample
        root_in.touch()
        root = root_in / f"experimental_data {sample}_{wavelength}.json"
        with open(root, "w") as w:
            json.dump(dictionary, w, cls=NumpyEncoder)

        logging.info("Data_proc:save_exp => ok !")

        return root

    def compil_experimental_data(self, sample, wavelength_list, root):
        """
        Unused function
        """

        logging.info("Data_proc:Starting compil_experimental_data.")
        logging.info(f"Data_proc:sample: {sample}.")
        logging.info(f"Data_proc:wavelength_list: {wavelength_list}.")
        logging.info(f"Data_proc:root: {root}.")

        root = root / sample
        temp_dict = {}
        for wavelength in wavelength_list:

            logging.info(f"Data_proc:wavelength: {wavelength}, wavelength_list: {wavelength_list}.")

            root_in = root / f"experimental_data {sample}_{wavelength}.json"
            with open(root_in, "r") as r:
                data = json.load(r)
            temp_dict.update({f"{wavelength}": data[f"{wavelength}"]})
            #os.remove(root_in)
        file = root / f"experimental_data {sample}.json"
        with open(file, "w") as w:
            json.dump(temp_dict, w, cls=NumpyEncoder)

        logging.info("Data_proc:compil_experimental_data => ok !")

        return file

    def save_process(self, dictionary, sample, state, simulation, root):
        """
        This function saves data in processed_data sample (temp).json file
        :param dictionary: dict
        :param sample: str
        :param state: bool
        :param simulation: bool
        :param root: path
        :return: path
        """

        logging.info("Data_proc:Starting save_process.")
        logging.info(f"Data_proc:type(dictionary): {type(dictionary)}.")
        logging.info(f"Data_proc:sample: {sample}.")
        logging.info(f"Data_proc:state: {state}.")
        logging.info(f"Data_proc:simulation: {simulation}.")
        logging.info(f"Data_proc:root: {root}.")

        if simulation == False:
            root = root/ f"{sample}" / f"processed_data {sample}.json"
        elif simulation == True:
            root = root / f"{sample}" / f"processed_data {sample}_temp.json"

        logging.info(f"Data_proc:simulation: {simulation}, state: {state},root: {root}.")

        if state == False :
            with open(root, "r") as r:
                old_data = json.load(r)
                old_data.update(dictionary)
            with open(root, "w") as w:
                json.dump(old_data, w, cls=NumpyEncoder)
        else:
            with open(root, "w") as w:
                json.dump(dictionary, w, cls=NumpyEncoder)

        logging.info("Data_proc:save_process => ok !")

        return root

    def read_process(self, sample, simulation, root):
        """
        This function reads data from processed_data sample (temp).json file
        :param sample: str
        :param simulation: bool
        :param root: path
        :return: dict
        """

        logging.info("Data_proc:Starting read_process.")
        logging.info(f"Data_proc:sample: {sample}.")
        logging.info(f"Data_proc:simulation: {simulation}.")
        logging.info(f"Data_proc:root: {root}.")

        if simulation == False:
            root = root / f"{sample}" / f"processed_data {sample}.json"
        elif simulation == True:
            root = root / f"{sample}" / f"processed_data {sample}_temp.json"
        with open(root, "r") as r:
            data = json.load(r)

        logging.info("Data_proc:read_process => ok !")

        return data

    def read_json_data(self, file):
        """
        This function retunr the data of the json input file.
        :param file: path
        :return: dict
        """

        logging.info(f"Data_proc:Starting read_json_data:file: {file}.")

        with open(file, "r") as r:
            data = json.load(r)

        logging.info("Data_proc:read_json_data => ok !")

        return data

    def save_sample(self, sample_name, emission_wavelength, emission_intensity, concentration, solvent, phi, root):
        """
        Unused function
        """

        logging.info("Data_proc:Starting save_sample.")
        logging.info(f"Data_proc:sample_name: {sample_name}.")
        logging.info(f"Data_proc:emission_wavelength: {emission_wavelength}.")
        logging.info(f"Data_proc:emission_intensity: {emission_intensity}.")
        logging.info(f"Data_proc:concentration: {concentration}.")
        logging.info(f"Data_proc:solvent: {solvent}.")
        logging.info(f"Data_proc:phi: {phi}.")
        logging.info(f"Data_proc:root: {root}.")

        data = {
            'solvent': pd.Series(solvent),
            'Concentration / M': pd.Series(concentration),
            'Phi': pd.Series(phi),
            'Reference wavelength / nm': pd.Series(emission_wavelength),
            'Reference intensity / cps': pd.Series(emission_intensity)
                }
        data_to_save = pd.DataFrame(data)
        root = root / f"{sample_name} informations.xlsx"
        with pd.ExcelWriter(root) as writer:
            data_to_save.to_excel(writer)

        logging.info("Data_proc:save_sample => ok !")

    def figure(self, x, y, name, sample, root):
        """
        This function draws figure containing a plot and save it as name.png to DATA/Figures folder.
        :param x: array
        :param y: array
        :param name: str
        :param sample: str
        :param root: path
        :return: bool
        """

        logging.info("Data_proc:Starting figure.")
        logging.info(f"Data_proc:type(x): {type(x)}.")
        logging.info(f"Data_proc:type(y): {type(y)}.")
        logging.info(f"Data_proc:name: {name}.")
        logging.info(f"Data_proc:sample: {sample}.")
        logging.info(f"Data_proc:root: {root}.")

        fig, ax = plt.subplots()
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter(
            '%.1e'))  # here %.1e is used to put scientific notation  with 1 significative value
        ax.xaxis.set_major_formatter(
            ticker.FormatStrFormatter('%.0f'))  # here %.0f is used to put decimal with 0 significative value
        ax.xaxis.set_major_locator(ticker.MaxNLocator(10))
        ax.xaxis.set_minor_locator(ticker.MaxNLocator(100))
        ax.yaxis.set_major_locator(ticker.MaxNLocator(6))
        ax.plot(x, y, color='black')
        plt.xlabel('Wavelength / nm')
        plt.ylabel('Fluorescence intensity / cps')
        DATA = root / "Figures"
        DATA.mkdir(exist_ok=True)
        DATA = root / "Figures" / "Experimental data"
        DATA.mkdir(exist_ok=True)
        DATA = root / "Figures" / "Experimental data" / f"{sample}"
        DATA.mkdir(exist_ok=True)
        fig.savefig(DATA / f"{name}.png")

        logging.info("Data_proc:figure => ok !")

        return True

    def figure_quadraticity(self, x1, y1, x2, y2, sample, wavelength, root):
        """
        This function draws figure containing a plot and save it as name.png to DATA/Figures folder.
        :param x1: array
        :param y1: array
        :param x2: array
        :param y2: array
        :param sample: str
        :param wavelength: int
        :param root: path
        :return: bool
        """

        logging.info("Data_proc:Starting figure_quadraticity.")
        logging.info(f"Data_proc:type(x1): {type(x1)}.")
        logging.info(f"Data_proc:type(y1): {type(y1)}.")
        logging.info(f"Data_proc:type(x2): {type(x2)}.")
        logging.info(f"Data_proc:type(y2): {type(y2)}.")
        logging.info(f"Data_proc:sample: {sample}.")
        logging.info(f"Data_proc:wavelength: {wavelength}.")
        logging.info(f"Data_proc:root: {root}.")

        fig, ax=plt.subplots()
        ax.plot(x1, y1, 'o', color='black', label='Original data')
        ax.plot(x2,y2, color='red', label='Fitted line')
        plt.xlabel('Wavelength / nm')
        plt.ylabel(r'Fluorescence intensity / cps')
        plt.xlim([-30,0])
        plt.ylim([0, 20])
        plt.legend()
        root = root / "Figures" / "Quadraticity"
        root.mkdir(exist_ok=True)
        root = root / f"{sample}"
        root.mkdir(exist_ok=True)
        fig.savefig(root / f"Quadraticity of {sample} at {wavelength} nm.png")

        logging.info("Data_proc:figure_quadraticity => ok !")

        return True

    def figure_emission_2p_vs_1p(self, x1, y1, x2, y2, sample, wavelength, root, power, wavenumber):
        """
        Unused function
        """

        logging.info("Data_proc:Starting figure_emission_2p_vs_1p.")
        logging.info(f"Data_proc:x1: {x1}.")
        logging.info(f"Data_proc:y1: {y1}.")
        logging.info(f"Data_proc:x2: {x2}.")
        logging.info(f"Data_proc:y2: {y2}.")
        logging.info(f"Data_proc:sample: {sample}.")
        logging.info(f"Data_proc:wavelength: {wavelength}.")
        logging.info(f"Data_proc:root: {root}.")
        logging.info(f"Data_proc:power: {power}.")
        logging.info(f"Data_proc:wavenumber: {wavenumber}.")

        fig, ax=plt.subplots()
        ax.plot(x1, y1, color='black', label='2P emission')
        ax.plot(x2, y2, color='red', label='1P emission')
        if wavenumber == False:
            plt.xlabel('Wavelength / nm')
        elif wavenumber == True:
            plt.xlabel('Wavenumber / cm-1')
        plt.ylabel(r'Fluorescence intensity / cps')
        #plt.xlim([-30,0])
        #plt.ylim([0, 20])
        plt.legend()
        root = root / "Figures" / "Experimental data" / f"{sample}" / f"{wavelength}"
        root.mkdir(exist_ok=True)
        fig.savefig(root / f"2P vs 1P emission of {sample} at {wavelength} nm _{power}.png")

        logging.info("Data_proc:figure_emission_2p_vs_1p => ok !")

        return True


    def save_experimental_to_excel(self, file, sample, root):
        """
        This function saves the experimental data of the sample in the corresponding xls file.
        :param file: path
        :param sample: str
        :param root: path
        :return: None
        """

        logging.info("Data_proc:Starting save_experimental_to_excel.")
        logging.info(f"Data_proc:file: {file}.")
        logging.info(f"Data_proc:sample: {sample}.")
        logging.info(f"Data_proc:root: {root}.")

        data = self.read_json_data(file)
        root = root / f"Experimental data of {sample}.xlsx"
        dictionary = {}
        test={}
        d = {}
        df = pd.DataFrame(d)
        for wavelength in data:
            dark = data[f"{wavelength}"]["dark"]
            for measure in data[f"{wavelength}"]:
                if measure != "dark":
                    emission_wavelength= data[f"{wavelength}"][f"{measure}"]["wavelength"]
                    intensity=data[f"{wavelength}"][f"{measure}"]["raw intensity"]
                    power=data[f"{wavelength}"][f"{measure}"]["raw power"]
                    intensity_corr= data[f"{wavelength}"][f"{measure}"]["intensity - dark"]
                    power_corr= data[f"{wavelength}"][f"{measure}"]["raw power - dark"]
                    area_corr= data[f"{wavelength}"][f"{measure}"]["fully corrected area"]
                    df["Dark - cps"] = pd.Series(dark)
                    df[f"{measure} - emission wavelength / nm"] = pd.Series(emission_wavelength)
                    df[f"{measure} - fluorescence intensity / cps"] = pd.Series(intensity)
                    df[f"{measure} - fluorescence intensity with dark correction / cps"] = pd.Series(intensity_corr)
                    df[f"{measure} - power"] = pd.Series(power)
                    df[f"{measure} - power with dark correction"] = pd.Series(power_corr)
                    df[f"{measure} - area with correction"] = pd.Series(area_corr)
                    dictionary.update({f"{measure}": df})
            temp = copy.deepcopy(dictionary)
            test[wavelength] = temp
        with pd.ExcelWriter(root) as writer:
            for key in test :
                for keys, value in test[key].items():
                    value.to_excel(writer, sheet_name=key)

        logging.info("Data_proc:save_experimental_to_excel => ok !")

    def save_process_data_to_excel(self, file, sample, root):
        """
        This function saves the processed data of the sample in the corresponding xls file.
        :param file: path
        :param sample: str
        :param root: path
        :return: None
        """

        logging.info("Data_proc:Starting save_process_data_to_excel.")
        logging.info(f"Data_proc:file: {file}.")
        logging.info(f"Data_proc:sample: {sample}.")
        logging.info(f"Data_proc:root: {root}.")

        data = self.read_json_data(file)
        dictionary = {}
        list_of_wavelength = []
        list_full_corr_area = []
        list_slope = []
        list_coeff= []
        list_log_fluorescence = []
        list_log_square_power=[]
        list_fit=[]
        list_S2F=[]
        list_S2=[]
        for wavelength in data:
            list_of_wavelength.append(wavelength)
            list_full_corr_area.append(data[f"{wavelength}"]["full_corrected_area"])
            list_slope.append(data[f"{wavelength}"]["quadraticity"]["slope"])
            list_coeff.append(data[f"{wavelength}"]["quadraticity"]["coeff"])
            list_log_fluorescence.append(data[f"{wavelength}"]["quadraticity"]["log(F)"])
            list_log_square_power.append(data[f"{wavelength}"]["quadraticity"]["log(P²)"])
            list_fit.append(data[f"{wavelength}"]["quadraticity"]["fit"])
            list_S2F.append(data[f"{wavelength}"]["S2F"])
            list_S2.append(data[f"{wavelength}"]["S2"])
        d = {
            "Emission wavelength / nm": list_of_wavelength,
            "S2F / GM": list_S2F,
            "S2 / GM": list_S2,
            "log(F)/log(P²)": list_slope,
            "fit of log(F)/log(P²)": list_fit,
            "R²":list_coeff,
            "log(F)": list_log_fluorescence,
            "log(P²)": list_log_square_power,
            "Full corrected area": list_full_corr_area
            }
        df = pd.DataFrame(d)
        root = root / f"Processed data of {sample}.xlsx"
        with pd.ExcelWriter(root) as writer:
            df.to_excel(writer)

        logging.info("Data_proc:save_process_data_to_excel => ok !")

    def figure_intensity_vs_power(self,number_of_measure, x, y, sample, wavelength, root):
        """
        This function draws figure containing a plot and save it.
        :param number_of_measure: int
        :param x: array
        :param y: array
        :param sample: str
        :param wavelength: int
        :param root: path
        :return: bool
        """

        logging.info("Data_proc:Starting figure_intensity_vs_power.")
        logging.info(f"Data_proc:number_of_measure: {number_of_measure}.")
        logging.info(f"Data_proc:type(x): {type(x)}.")
        logging.info(f"Data_proc:type(y): {type(y)}.")
        logging.info(f"Data_proc:sample: {sample}.")
        logging.info(f"Data_proc:wavelength: {wavelength}.")
        logging.info(f"Data_proc:root: {root}.")

        dictionary_color = {"0": "C0", "1": "C1", "2": "C2", "3": "C3", "4": "C4", "5": "C5", "6": "C6", "7": "C7", "8": "C8", "9": "C9"}
        fig, ax = plt.subplots()
        for measure in range(number_of_measure+1):
            if measure !=0:
                if measure > 9:
                    data = str(measure)
                    data_temp = data[1:]
                    measure = int(data_temp)
                ax.plot(x[measure], y[measure], color=dictionary_color[f"{measure}"])
            else:
                pass
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter(
            '%.1e'))  # here %.1e is used to put scientific notation  with 1 significative value
        ax.xaxis.set_major_formatter(
            ticker.FormatStrFormatter('%.0f'))  # here %.0f is used to put decimal with 0 significative value
        ax.xaxis.set_major_locator(ticker.MaxNLocator(10))
        ax.xaxis.set_minor_locator(ticker.MaxNLocator(100))
        ax.yaxis.set_major_locator(ticker.MaxNLocator(6))
        plt.xlabel('Wavelength / nm')
        plt.ylabel('Fluorescence intensity / cps')
        root = root / "Figures" / "Experimental data" / f"{sample}"
        root.mkdir(exist_ok=True)
        fig.savefig(root / f"{sample} at {wavelength} nm.png")

        logging.info("Data_proc:figure_intensity_vs_power => ok !")
        return True
