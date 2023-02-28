"""
Created on 5/1/2023
@author: Jonathan DANIEL, jonathan.daniel@u-bordeaux.fr
"""
SIMULATION_FOR_DEBUG= False  #Set True only for code test and debug, for real measure set False
                            # do not forget to do the same for PyTPEF

if SIMULATION_FOR_DEBUG== True:
    import Ressources_scripts.Spectro_Simulation as talk_to_spectro
    #import Ressources_scripts.Photodiode_Simulation as nidaqmx
else:
    import Ressources_scripts.talk_to_spectro as talk_to_spectro
    #import nidaqmx
    #from Ressources_scripts import PMA100
    import Ressources_scripts.Photodiode as Photodiode


import copy
import os
import json
import time
from pathlib import Path
import pandas as pd
from threading import Thread
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
temp = Path(f"{CUR_DIR}").parent
temp = temp / "Ressources" / "save_directory.json"

with open(temp, "r") as r:
    save_dir= json.load(r)
    print(save_dir)
DATA = Path(f"{save_dir}/Datas")
DATA.mkdir(exist_ok=True)
exp = DATA / "Figures"
exp.mkdir(exist_ok=True)
exp_2=exp / "Experimental data"
exp_2.mkdir(exist_ok=True)
ROOT = DATA / "Metadata"
ROOT.mkdir(exist_ok=True)
root= ROOT / "sample_data.json"
root.touch()



def get_path():
    """
    This function only return Path of Data folder
    :return: Path
    """

    logging.info(f"Config:Starting get_path: DATA: {DATA}.")

    return DATA

def create_json():
    """
    This function create the file sample_data.json
    :return: None
    """

    logging.info(f"Config:Starting create_json.")

    r = "[]"
    with open(root, "w") as g:
        json.dump(r, g)

    logging.info("Config:create_json => ok !")

def create_samples(list_of_samples):
    """
    This function saves the list of the samples name in DATA/samples.json file
    :param list_of_samples: lst
    :return:None
    """

    logging.info(f"Config:Starting create_samples.")
    logging.info(f"Config:list_of_samples: {list_of_samples}.")

    root_in = ROOT / "samples.json"
    with open(root_in, "w") as w:
        json.dump(list_of_samples, w)

    logging.info("Config:create_samples => ok !")

def get_samples():
    """
    This function just return the list of all samples store in DATA/samples.json file
    :return: lst
    """
    logging.info(f"Config:Starting get_samples.")

    samples= []
    root_in = ROOT / "samples.json"
    with open(root_in, "r") as f:
        sample_names = json.load(f)
        for _ in sample_names:
            samples.append(sample_data(_))

        logging.info("Config:get_samples => ok !")

        return samples

def get_wavelength(sample):
    """
    Unused function
    """

    logging.info(f"Config:Starting get_wavelength.")
    logging.info(f"Config:sample: {sample}.")

    wavelength = []
    root_in = DATA / f"{sample}" / "wavelength.json"
    with open(root_in, "r") as f:
        wavelengths = json.load(f)
        for _ in wavelengths:
            wavelength.append(sample_data(str(_)))

        logging.info("Config:get_wavelength => ok !")

        return wavelength

class Temp_Files():

    def __init__(self):
        pass

    def save_files(self, data, name):
        """
        Unused function
        """

        logging.info(f"Config:Temp_Files:Starting save_files.")
        logging.info(f"Config:Temp_Files:data: {data}.")
        logging.info(f"Config:Temp_Files:name: {name}.")

        root_in = DATA / "temp"
        root_in.mkdir(exist_ok=True)
        root_in = root_in / f"{name}.json"
        root_in.touch()
        if name =="root_dict":
            data= f"{data}"
        with open(root_in, "w") as w:
            json.dump(data, w, cls= NumpyEncoder)

        logging.info("Config:Temp_Files:save_files => ok !")

    def get_files(self, name):
        """
        Unused function
        """

        logging.info(f"Config:Temp_Files:Starting get_files.")
        logging.info(f"Config:Temp_Files:name: {name}.")

        root_in = DATA / "temp" / f"{name}.json"
        with open(root_in, "r") as r:
            data = json.load(r)

        logging.info("Config:Temp_Files:get_files => ok !")

        return data

class chrono:

    def __init__(self):
        pass

    def run(self):

        logging.info(f"Config:chrono:Starting run.")

        self.start = time.time()

        logging.info("Config:chrono:run => ok !")

        return self.start

    def restart(self):

        logging.info(f"Config:chrono:Starting restart.")

        self.start = time.time()

        logging.info("Config:chrono:restart => ok !")

    def get_new_time(self, start):
        """
        Unused function
        """

        logging.info(f"Config:chrono:Starting get_new_time.")
        logging.info(f"Config:chrono:start: {start}.")

        value = time.time() - start
        self.restart()

        logging.info("Config:chrono:get_new_time => ok !")

        return value

    def get_time(self, start):

        logging.info(f"Config:chrono:Starting get_time.")
        logging.info(f"Config:chrono:start: {start}.")

        actual_time = time.time()-start

        logging.info("Config:chrono:get_time => ok !")

        return actual_time

class parrallel_execution:

    def __init__(self, integration_time):

        # try:
        #     # self.power_meter = PMA100.PMA100()
        #     # self.power_meter.Connect()
        #     # logging.info("Connected to PMA100")
        # except:
        #     logging.error("CANNOT CONNECT TO PMA100 DEVICE !")

        logging.info(f"Config:parrallel_execution:integration_time: {integration_time}.")

        self.integration_time = integration_time
        self.power_data=[]
        self.spectro= talk_to_spectro.Maya(self.integration_time)
        self.repeat = True

        logging.info("Config:parrallel_execution:Initialization => ok !")

    def _spectro(self):

        logging.info(f"Config:parrallel_execution:Starting _spectro.")
        logging.info("Config:Start acquisition")

        self.wavelength, self.spectrum, self.raw_fluorescence=self.spectro._data_read_scan(self.scan, self.intensity_dark)
        self.repeat = False

        logging.info("Config:Step 1 - Spectra acquired")

        time.sleep(0.5)

        logging.info("Config:End of acquisition")
        logging.info(f"Config:Length of power_data: {len(self.power_data)}, Integration time: {self.integration_time}.")
        logging.info("Config:parrallel_execution:_spectro => ok !")

    def _photodiode(self):

        logging.info(f"Config:parrallel_execution:Starting _photodiode.")
        logging.info("Config:Please Remove the comments on _photodiode, line 136 for real hardware !!")

        while self.repeat == True :
            if SIMULATION_FOR_DEBUG :
                data=np.random.randint(1,10)
                self.power_data.append(data)
            else:
                # time, data = self.power_meter.Read_Power()
                data = Photodiode.read()*1000 # W to mW
                """with nidaqmx.Task() as task:
                    task.ai_channels.add_ai_voltage_chan("Dev1/ai1")
                    data_power = task.read()
                    self.power_data.append(data_power)
                task.close()"""
                if str(data) == "-inf" or str(data) == "inf":
                    pass
                else:
                    self.power_data.append(data)

        logging.info("Config:Power measure ended.")
        logging.info("Config:parrallel_execution:_photodiode => ok !")

    def parallel_run(self):

        logging.info(f"Config:parrallel_execution:Starting parallel_run.")

        a = Thread(target=self._spectro)
        b = Thread(target=self._photodiode)
        a.start()
        b.start()
        """a.join() is thread function which block the script while the thread is running"""
        a.join()

        logging.info("Config:Thread.join() successful !")
        logging.info("Config:parrallel_execution:parallel_run => ok !")

    def run(self, scan, intensity_dark):

        logging.info(f"Config:parrallel_execution:Starting run.")
        logging.info(f"Config:parrallel_execution:scan: {scan}.")
        logging.info(f"Config:parrallel_execution:intensity_dark: {intensity_dark}.")

        self.scan= scan
        self.intensity_dark = intensity_dark
        self.parallel_run()

        logging.info("Config:parrallel_execution:run => ok !")
        # self.power_meter.Disconnect()
        return self.wavelength, self.spectrum, self.power_data, self.raw_fluorescence

class sample_data:

    def __init__(self, sample_name):
        """
        the file DATA/sample_data.json is a buffer file for storing samples names & data.
        :param sample_name:str
        """

        logging.info(f"Config:sample_data:sample_name: {sample_name}.")

        self.sample_name = sample_name
        with open(root, 'w') as f:
            json.dump([], f)

        logging.info("Config:sample_data:Initialization => ok !")

    def __str__(self):
        return self.sample_name

    def _get_sample_data(self, remove=False):
        """
        This function return the data stored in the buffer file.
        :return: str
        """

        logging.info(f"Config:sample_data:Starting _get_sample_data.")
        if remove:
            root = ROOT / "samples.json"
        else:
            root = ROOT / "sample_data.json"
        with open(root, "r") as f:

            logging.info("Config:sample_data:_get_sample_data => ok !")

            return json.load(f)

    def _write_sample_data(self, sample_name, remove=False):
        """
        This function write sample_name in the buffer file.
        :param sample_name: str
        :return: None
        """

        logging.info(f"Config:sample_data:Starting _write_sample_data.")
        logging.info(f"Config:sample_data:sample_name: {sample_name}.")


        if remove:
            root = ROOT / "samples.json"
        else:
            root = ROOT / "sample_data.json"
        with open(root, 'w') as f:
            json.dump(sample_name, f)

            logging.info("Config:sample_data:_write_sample_data => ok !")

    def add_to_sample_name(self):
        """
        :return:bool
        """

        logging.info(f"Config:sample_data:Starting add_to_sample_name.")

        data = self._get_sample_data()
        data.append(self.sample_name)

        logging.info(f"Config:sample_data:data: {data}.")

        self._write_sample_data(data)

        logging.info("Config:sample_data:add_to_sample_name => ok !")

        return True

    def remove_sample(self):
        """
        This function remove the sample from buffer file
        :return:
        """

        logging.info(f"Config:sample_data:Starting remove_sample {self.sample_name}.")

        data = self._get_sample_data(remove=True)

        logging.info(f"Config:sample_data:data: {data}.")

        if self.sample_name in data:
            data.remove(self.sample_name)
            self._write_sample_data(data, remove=True)

        logging.info("Config:sample_data:remove_sample => ok !")

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
        self.root = DATA
        self.dict_angle_power_root = self.root / "Figures" / "Experimental data" / "data_power_angle.json"

    def create_angle_power_dict_for_measure(self):
        """
        This function create a file (see init) used to assign angle of half-wave plate to real laser power of angle
        used for the 2PA measurement.
        :return:None
        """

        logging.info(f"Config:save:Starting create_angle_power_dict_for_measure.")

        self.dict_angle_power_root.touch()

        logging.info("Config:save:create_angle_power_dict_for_measure => ok !")

    def fill_angle_power_dict_for_measure(self, data):
        """
        This function write the data inside the angle-power file.
        :param data:dict
        :return: bool
        """

        logging.info(f"Config:save:Starting fill_angle_power_dict_for_measure.")
        logging.info(f"Config:save:data: {data}.")

        with open(self.dict_angle_power_root, "w") as w:
            json.dump(data, w)

        logging.info("Config:save:fill_angle_power_dict_for_measure => ok !")

        return True

    def read_angle_power_dict_for_measure(self):
        """
        This function read the data inside the angle-power file.
        :return: dict
        """

        logging.info(f"Config:save:Starting read_angle_power_dict_for_measure.")

        with open(self.dict_angle_power_root, "r") as r:
            data= json.load(r)

        logging.info("Config:save:read_angle_power_dict_for_measure => ok !")

        return data

    def extract_angle_power_from_dict(self, wavelength, i):
        """
                This function extract data from the angle-power file.
        :param wavelength: int
        :param i: int (position of the power measure)
        :return: float, float
        """

        logging.info(f"Config:save:Starting extract_angle_power_from_dict.")
        logging.info(f"Config:save:wavelength: {wavelength}.")
        logging.info(f"Config:save:i: {i}.")

        data=self.read_angle_power_dict_for_measure()
        get_data=data[f"{wavelength}"][f"{i}"]
        angle=get_data[0]
        power=get_data[1]

        logging.info(f"Config:Angle to measure again is: {angle}, Power to measure again is {power}.")
        logging.info("Config:save:extract_angle_power_from_dict => ok !")

        return angle, power


    def fit_curve(self, x, y):
        """
        Unused function
        """
        """
        This function fits the curve log(y) =f(log(x²)) using a 2nd order polynome.
        The slope corresponds to the order of the absorption process
        :param x: array
        :param y: array
        :return: array, float, float
        """

        logging.info(f"Config:save:Starting fit_curve.")
        logging.info(f"Config:save:type(x): {type(x)}.")
        logging.info(f"Config:save:type(y): {type(y)}.")

        def fit_function(x, a, b):

            logging.info(f"Config:save:fit_function.")
            logging.info(f"Config:save:type(x): {type(x)}.")
            logging.info(f"Config:save:a: {a}.")
            logging.info(f"Config:save:b: {b}.")

            return ((a * x) + b)

        try:

            # curve_fit() function takes the test-function
            # x-data and y-data as argument and returns
            # the coefficients a and b in param and
            # the estimated covariance of param in param_cov
            param, param_cov = curve_fit(fit_function, x, y, p0=[0.5,0.5])
            line= np.polyfit(x, y, 1)
            #line = np.poly1d(np.polyfit(x, y, 1))

            logging.info(f"Config:line: {line}.")

            lines=line[0]*x+line[1]

            logging.info("Config:function coefficients:")
            logging.info(f"Config:param: {param}.")
            logging.info("Config:Covariance of coefficients:")
            logging.info(f"Config:param_cov: {param_cov}.")

            # ans stores the new y-data according to
            # the coefficients given by curve-fit() function

            logging.info(f"Config:param[0]: {param[0]}.")
            logging.info(f"Config:param[1]: {param[1]}.")

            ans = ((param[0] * x) + param[1])

            """plt.plot(x, y, 'o', color ='red', label ="data")
            plt.plot(x, lines, '--', color='green', label="optimized data with line")
            plt.legend()
            plt.show()"""

            reg_lin= linregress(np.log(x), np.log(y))
            slope= round(reg_lin[0],1)
            correlation_coeff=round(reg_lin[2]*reg_lin[2],5)

        except RuntimeError:
            lines=2*x + 0
            correlation_coeff = 0
            slope = 0

            logging.error("Config:Error, cannot fit the data !!")
        logging.info("Config:save:fit_curve => ok !")

        return lines, correlation_coeff, slope

    def create_samples_informations(self, data_to_save):
        """
        This function save metadata concerning the samples in the file samples and measure informations.json.
        :param data_to_save: dict
        :return: None
        """

        logging.info(f"Config:save:Starting create_samples_informations.")
        logging.info(f"Config:save:data_to_save: {data_to_save}.")

        root_in =self.root / "Figures" / "Experimental data" / "samples and measure informations.json"
        root_in.touch()
        with open(root_in, "w") as w:
            json.dump(data_to_save, w)

        logging.info("Config:save:create_samples_informations => ok !")

    def extract_process_data(self, sample, wavelength):
        """
        Unused function
        """
        """
        This function read the processed data of the sample at specific wavelength.
        :param sample: str
        :param wavelength: int
        :return:lst, lst, lst, lst
        """

        logging.info(f"Config:save:Starting extract_process_data.")
        logging.info(f"Config:save:sample: {sample}.")
        logging.info(f"Config:save:wavelength: {wavelength}.")

        list_of_s2f = []
        list_of_wavelength = []
        root_in = DATA / f"{sample}" / f"processed_data {sample}.json"
        with open(root_in, "r") as r:
            data = json.load(r)

            logging.info(f"Config:extract_process_data: {data.keys()}.")

        for key_wavelength in data:
            temp_var=data[f"{key_wavelength}"].keys()

            logging.info(f"Config:extract_process_data_bis: {temp_var}.")

            list_of_s2f.append(data[f"{key_wavelength}"]["S2F"])
            list_of_wavelength.append(key_wavelength)
        if wavelength == None:
            list_of_log_f = False
            list_of_log_pp = False
        else:
            list_of_log_f= data[f"{wavelength}"]["quadraticity"]["log(F)"]
            list_of_log_pp= data[f"{wavelength}"]["quadraticity"]["log(P\u00b2)"]

        logging.info("Config:save:extract_process_data => ok !")

        return list_of_s2f, list_of_wavelength, list_of_log_f, list_of_log_pp

    def remove_data(self, sample, wavelength, measure):
        """
        Unused function
        """
        """
        This function remove the data of a specific measure at a given wavelength for the sample.
        :param sample: str
        :param wavelength: int
        :param measure: int
        :return: None
        """

        logging.info(f"Config:save:Starting remove_data.")
        logging.info(f"Config:save:sample: {sample}.")
        logging.info(f"Config:save:wavelength: {wavelength}.")
        logging.info(f"Config:save:measure: {measure}.")

        root_in = DATA / f"{sample}" / f"experimental_data {sample}_{wavelength}.json"
        with open(root_in, "r") as r:
            data = json.load(r)
        data[f"{wavelength}"].pop(f"{measure}")

        logging.info("Config:Measure delete successfuly !")

        with open(root_in, "w") as w:
            json.dump(data, w)

        logging.info("Config:save:remove_data => ok !")

    def _get_file_for_new_mesure(self, sample, wavelength):
        """
        This function give the path for a new measure of the sample at a given wavelength.
        :param sample: str
        :param wavelength: int
        :return: path
        """

        logging.info(f"Config:save:Starting _get_file_for_new_mesure.")
        logging.info(f"Config:save:sample: {sample}.")
        logging.info(f"Config:save:wavelength: {wavelength}.")

        root_in = DATA / f"{sample}" / f"experimental_data {sample}_{wavelength}.json"

        logging.info("Config:save:_get_file_for_new_mesure => ok !")

        return root_in

    def _get_data_for_new_mesure(self, sample, wavelength):
        """
        Unused function
        """
        """
        This function gives the data of a new measure of the sample at a given wavelength.
        :param sample: str
        :param wavelength: int
        :return: dict
        """

        logging.info(f"Config:save:Starting _get_data_for_new_mesure.")
        logging.info(f"Config:save:sample: {sample}.")
        logging.info(f"Config:save:wavelength: {wavelength}.")

        root = self._get_file_for_new_mesure(sample= sample, wavelength= wavelength)
        with open(root, "r") as r:
            data = json.load(r)

        logging.info("Config:save:_get_data_for_new_mesure => ok !")

        return data

    def create_file(self, sample):
        """
        This function creates the file DATA/sample/wavelength.json and fill it with an empty list
        :param sample: str
        :return: None
        """

        logging.info(f"Config:save:Starting create_file.")
        logging.info(f"Config:save:sample: {sample}.")

        first_root = DATA / f"{sample}"
        first_root.mkdir(exist_ok=True)
        root_in = DATA / f"{sample}" / "wavelength.json"
        root_in.touch()
        void = []
        with open(root_in, "w") as w:
            json.dump(void, w)

        logging.info("Config:save:create_file => ok !")

    def _get_location_file(self):
        """
        Unused function
        """
        """
        This function return the path DATA
        :return: path
        """

        logging.info(f"Config:save:_get_location_file.")

        return DATA

    def write_only(self, dictionary_process, sample):
        """
        This function write the data in the processed_data sample.json file.
        :param dictionary_process: dict
        :param sample: str
        :return: None
        """

        logging.info(f"Config:save:Starting write_only.")
        logging.info(f"Config:save:type(dictionary_process): {type(dictionary_process)}.")
        logging.info(f"Config:save:sample: {sample}.")

        root_in = DATA / sample
        root_in.mkdir(exist_ok=True)
        root = root_in / f"processed_data {sample}.json"
        root.touch()
        with open(root, "w") as w:
            json.dump(dictionary_process, w, cls=NumpyEncoder)

        logging.info("Config:save:write_only => ok !")

    def extract_datas(self, dictionnary):
        """
        Unused function
        """
        """
        This function extract data from the input dictionnary
        :param dictionnary: dict
        :return: array, array
        """

        logging.info(f"Config:save:Starting extract_datas.")
        logging.info(f"Config:save:type(dictionnary): {type(dictionnary)}.")

        x = []
        y = []
        for key in dictionnary:
            value = dictionnary[key]["quadraticity"]["slope"]
            x.append(key)
            y.append(value)
        x = np.array(x)
        y = np.array(y)

        logging.info("Config:save:extract_datas => ok !")

        return x, y

    def save_exp(self, dictionary, sample, wavelength):
        """
        This function save the data in experimental_data sample_wavelength.json file
        :param dictionary: dict
        :param sample: str
        :param wavelength: int
        :return: path
        """

        logging.info(f"Config:save:Starting save_exp.")
        logging.info(f"Config:save:type(dictionary): {type(dictionary)}.")
        logging.info(f"Config:save:sample: {sample}.")
        logging.info(f"Config:save:wavelength: {wavelength}.")

        root_in = DATA / f"{sample}"
        root_in.mkdir(exist_ok=True)
        root = root_in / f"experimental_data {sample}_{wavelength}.json"
        with open(root, "w") as w:
            json.dump(dictionary, w, cls=NumpyEncoder)

        logging.info("Config:save:save_exp => ok !")

        return root

    def compil_experimental_data(self, sample, wavelength_list):
        """
        Unused function
        """
        """
        This function read all the experimental data of the sample and put them in
        the experimental_data sample.json file
        :param sample: str
        :param wavelength_list: lst
        :return: path
        """

        logging.info(f"Config:save:Starting compil_experimental_data.")
        logging.info(f"Config:save:sample: {sample}.")
        logging.info(f"Config:save:wavelength_list: {wavelength_list}.")

        root = DATA / sample
        temp_dict = {}
        for wavelength in wavelength_list:

            logging.info(f"Config:wavelength: {wavelength}, wavelength_list: {wavelength_list}.")

            root_in = root / f"experimental_data {sample}_{wavelength}.json"
            with open(root_in, "r") as r:
                data = json.load(r)
            temp_dict.update({f"{wavelength}": data[f"{wavelength}"]})
            #os.remove(root_in)
        file = root / f"experimental_data {sample}.json"
        with open(file, "w") as w:
            json.dump(temp_dict, w, cls=NumpyEncoder)

        logging.info("Config:save:compil_experimental_data => ok !")

        return file

    def save_process(self, dictionary, sample, state):
        """
        This function save or update data (depends on state) of the sample in the processed_data sample.json file
        :param dictionary: dict
        :param sample: str
        :param state: bool
        :return: path
        """

        logging.info(f"Config:save:Starting save_process.")
        logging.info(f"Config:save:type(dictionary): {type(dictionary)}.")
        logging.info(f"Config:save:sample: {sample}.")
        logging.info(f"Config:save:state: {state}.")

        root = DATA/ f"{sample}" / f"processed_data {sample}.json"
        if state == False :
            with open(root, "r") as r:
                old_data = json.load(r)
                old_data.update(dictionary)
            with open(root, "w") as w:
                json.dump(old_data, w, cls=NumpyEncoder)
        else:
            with open(root, "w") as w:
                json.dump(dictionary, w, cls=NumpyEncoder)

        logging.info("Config:save:save_process => ok !")

        return root

    def read_process(self, sample):
        """
        This function read the processed_data sample.json file and return the data.
        :param sample: str
        :return: dict
        """

        logging.info(f"Config:save:Starting read_process.")
        logging.info(f"Config:save:sample: {sample}.")

        root = DATA / f"{sample}" / f"processed_data {sample}.json"
        with open(root, "r") as r:
            data = json.load(r)

        logging.info("Config:save:read_process => ok !")

        return data

    def read_json_data(self, file):
        """
        This function only read the json fil input and return the output data
        :param file: path
        :return: dict
        """

        logging.info(f"Config:save:Starting read_json_data.")
        logging.info(f"Config:save:file: {file}.")

        with open(file, "r") as r:
            data = json.load(r)

        logging.info("Config:save:read_json_data => ok !")

        return data

    def save_power_values(self, sample_name, excitation_wavelength, list_of_power_to_save):
        """
        This function saves exclusively all the power values of the sample at specific wavelength in a .txt file
        :param sample_name: str
        :param excitation_wavelength: int
        :param list_of_power_to_save: dict
        """

        file = self.root / f"{sample_name}" / f"{excitation_wavelength}_power_values.txt"
        file.touch()
        i=0
        head = "# \t"
        for key in list_of_power_to_save:
            head = head + f" \t {key}"
        head = head + "\n"
        with open(file, "w") as w:
            w.write(head)
        for n in list_of_power_to_save['dark']:
            try:
                data =f"{i} \t"
                for key in list_of_power_to_save:
                    value = list_of_power_to_save[f"{key}"][i]
                    data = data + f" \t {value}"
                data = data + "\n"
                with open(file, "a") as a:
                    a.write(data)
                i+=1
            except:
                pass

    def save_sample(self, sample_name, emission_wavelength, emission_intensity, concentration, solvent, phi):
        """
        This function saves the parameter data of the sample in the sample information.xls file
        :param sample_name: str
        :param emission_wavelength: array
        :param emission_intensity: array
        :param concentration: float
        :param solvent: str
        :param phi: float
        :return: None
        """

        logging.info(f"Config:save:Starting save_sample.")
        logging.info(f"Config:save:sample_name: {sample_name}.")
        logging.info(f"Config:save:emission_wavelength: {emission_wavelength}.")
        logging.info(f"Config:save:emission_intensity: {emission_intensity}.")
        logging.info(f"Config:save:concentration: {concentration}.")
        logging.info(f"Config:save:solvent: {solvent}.")
        logging.info(f"Config:save:phi: {phi}.")

        data = {
            'solvent': pd.Series(solvent),
            'Concentration / M': pd.Series(concentration),
            'Phi': pd.Series(phi),
            'Reference wavelength / nm': pd.Series(emission_wavelength),
            'Reference intensity / cps': pd.Series(emission_intensity)
                }
        data_to_save = pd.DataFrame(data)
        root = self.root / f"{sample_name} informations.xlsx"
        with pd.ExcelWriter(root) as writer:
            data_to_save.to_excel(writer)

        logging.info("Config:save:save_sample => ok !")

    def figure(self, x, y, name, sample):
        """
        This function draws figure containing a plot and save it as name.png to DATA/Figures folder.
        :param x: array
        :param y: array
        :param name: str
        :param sample: str
        :return: bool
        """

        logging.info(f"Config:save:Starting figure.")
        logging.info(f"Config:save:type(x): {type(x)}.")
        logging.info(f"Config:save:type(y): {type(y)}.")
        logging.info(f"Config:save:name: {name}.")
        logging.info(f"Config:save:sample: {sample}.")

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
        DATA = self.root / "Figures"
        DATA.mkdir(exist_ok=True)
        DATA = self.root / "Figures" / "Experimental data"
        DATA.mkdir(exist_ok=True)
        DATA = self.root / "Figures" / "Experimental data" / f"{sample}"
        DATA.mkdir(exist_ok=True)
        fig.savefig(DATA / f"{name}.png")

        logging.info("Config:save:figure => ok !")

        return True

    def figure_quadraticity(self, x1, y1, x2, y2, sample, wavelength):
        """
        This function draws figure containing a plot and save it as name.png to DATA/Figures folder.
        :param x1: array
        :param y1: array
        :param x2: array
        :param y2: array
        :param sample: str
        :param wavelength: int
        :return: bool
        """

        logging.info(f"Config:save:Starting figure_quadraticity.")
        logging.info(f"Config:save:type(x1): {type(x1)}.")
        logging.info(f"Config:save:type(y1): {type(y1)}.")
        logging.info(f"Config:save:type(x2): {type(x2)}.")
        logging.info(f"Config:save:type(y2): {type(y2)}.")
        logging.info(f"Config:save:sample: {sample}.")
        logging.info(f"Config:save:wavelength: {wavelength}.")

        fig, ax=plt.subplots()
        ax.plot(x1, y1, 'o', color='black', label='Original data')
        ax.plot(x2,y2, color='red', label='Fitted line')
        plt.xlabel('Wavelength / nm')
        plt.ylabel(r'Fluorescence intensity / cps')
        plt.xlim([-30,0])
        plt.ylim([0, 20])
        plt.legend()
        root = self.root / "Figures" / "Quadraticity"
        root.mkdir(exist_ok=True)
        root = self.root / "Figures" / "Quadraticity" / f"{sample}"
        root.mkdir(exist_ok=True)
        fig.savefig(root / f"Quadraticity of {sample} at {wavelength} nm.png")

        logging.info("Config:save:figure_quadraticity => ok !")

        return True

    def save_experimental_to_excel(self, file, sample):
        """
        Unused function
        """
        """
        This function saves the experimental data of the sample in the corresponding xls file.
        :param file: path
        :param sample: str
        :return: None
        """

        logging.info(f"Config:save:Starting create_angle_power_dict_for_measure.")
        logging.info(f"Config:save:file: {file}.")
        logging.info(f"Config:save:sample: {sample}.")

        data = self.read_json_data(file)
        root = self.root / f"Experimental data of {sample}.xlsx"
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

        logging.info("Config:save:create_angle_power_dict_for_measure => ok !")

    def save_process_data_to_excel(self, file, sample):
        """
        This function saves the processed data of the sample in the corresponding xls file.
        :param file: path
        :param sample: str
        :return: None
        """

        logging.info(f"Config:save:Starting save_process_data_to_excel.")
        logging.info(f"Config:save:file: {file}.")
        logging.info(f"Config:save:sample: {sample}.")

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
        root = self.root / f"Processed data of {sample}.xlsx"
        with pd.ExcelWriter(root) as writer:
            df.to_excel(writer)

        logging.info("Config:save:save_process_data_to_excel => ok !")

    def save_angle_to_power_conversion_to_excel(self, file):
        """
        Unused function
        """
        """
        This function saves the used angle and real power in Experimental measure of power vs angle of
        half-wave plate.xlsx file.
        :param file: path
        :return: None
        """

        logging.info(f"Config:save:Starting save_angle_to_power_conversion_to_excel.")
        logging.info(f"Config:save:file: {file}.")

        data = self.read_json_data(file)
        dictionary = {}
        power=[]
        angle=[]
        for wavelength in data:
            for key, values in data[f"{wavelength}"].items():
                power.append(values)
                angle.append(key)
            d = {
                f"{wavelength} - angle": angle,
                f"{wavelength} - power": power
            }
            df = pd.DataFrame(d)
            dictionary.update({f"{wavelength}":df})
            power.clear()
            angle.clear()
            logging.info(f"wavelength: {wavelength}.")
        root = self.root / f"Experimental measure of power vs angle of half-wave plate.xlsx"
        with pd.ExcelWriter(root) as writer:
            for key, value in dictionary.items():
                value.to_excel(writer, sheet_name=key)

        logging.info("Config:save:save_angle_to_power_conversion_to_excel => ok !")

    def figure_intensity_vs_power(self,number_of_measure, x, y, sample, wavelength):
        """
        This function draws figure containing a plot and save it.
        :param number_of_measure: int
        :param x: array
        :param y: array
        :param sample: str
        :param wavelength: int
        :return: bool
        """

        logging.info(f"Config:save:Starting figure_intensity_vs_power.")
        logging.info(f"Config:save:number_of_measure: {number_of_measure}.")
        logging.info(f"Config:save:type(x): {type(x)}.")
        logging.info(f"Config:save:type(y): {type(y)}.")
        logging.info(f"Config:sample:y2: {sample}.")
        logging.info(f"Config:save:wavelength: {wavelength}.")

        dictionary_color = {"0": "C0", "1": "C1", "2": "C2", "3": "C3", "4": "C4", "5": "C5", "6": "C6", "7": "C7", "8": "C8", "9": "C9"}
        fig, ax = plt.subplots()
        for measure in range(number_of_measure+1):
            if measure !=0:
                ax.plot(x[measure], y[measure], color=dictionary_color[f"{measure}"])
                logging.info(f"x[measure] is: {x[measure]}, type(x[measure]): {type(x[measure])}.")
                logging.info(f"y[measure] is: {y[measure]}, type(y[measure]): {type(y[measure])}.")
                logging.info(f"color is: {dictionary_color[measure]}.")
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
        root = self.root / "Figures" / "Experimental data" / f"{sample}"
        root.mkdir(exist_ok=True)
        fig.savefig(root / f"{sample} at {wavelength} nm.png")

        logging.info("Config:save:figure_intensity_vs_power => ok !")

        return True

class Config:

    def __init__(self):
        self.wavelength_range= []
        self.current_wavelength = None
        self.wavelength_min = None
        self.wavelength_max = None
        self.wavelength_pitch = None
        self.power_range= []
        self.power_min = None
        self.power_max = None
        self.power_pitch = None

    def refresh(self):
        return

    def wavelength(self, min, pitch, max):
        """
        This function creates the list of wavelength range to measure.
        :param min: int
        :param pitch: int
        :param max: int
        :return: lst
        """

        logging.info(f"Config:Starting wavelength.")
        logging.info(f"Config:min: {min}.")
        logging.info(f"Config:pitch: {pitch}.")
        logging.info(f"Config:max: {max}.")

        self.wavelength_range.clear
        self.wavelength_min = min
        self.wavelength_max = max
        self.wavelength_pitch = pitch
        self.wavelength_range.append(self.wavelength_min)
        self.current_wavelength = min
        if self.wavelength_max > self.wavelength_min:
            while self.current_wavelength < self.wavelength_max:
                self.current_wavelength += self.wavelength_pitch
                self.wavelength_range.append(self.current_wavelength)
        else:
            while self.current_wavelength > self.wavelength_max:
                self.current_wavelength -= self.wavelength_pitch
                self.wavelength_range.append(self.current_wavelength)

        logging.info(f"Config:self.wavelength_range: {self.wavelength_range}.")
        logging.info("Config:wavelength => ok !")

        return self.wavelength_range


    def power(self, min, number_of_measure, max):
        """
        Unused function
        """
        """
        This function creates the list of power range to measure.
        :param min: int
        :param number_of_measure: int
        :param max: int
        :return: lst
        """

        logging.info(f"Config:Starting power.")
        logging.info(f"Config:min: {min}.")
        logging.info(f"Config:number_of_measure: {number_of_measure}.")
        logging.info(f"Config:max: {max}.")

        self.power_min = min
        self.power_max = max
        self.power_number = number_of_measure
        self.current_power = min

        power = self.power_max - self.power_min
        power /= self.power_number
        self.power_range.append(self.power_min)

        if self.power_max > self.power_min :
            while self.current_power < self.power_max:
                self.current_power += power
                self.power_range.append(self.current_power)
        else:
            while self.current_power > self.power_max:
                self.current_power += power
                self.power_range.append(self.current_power)

        logging.info(f"Config:self.power_range: {self.power_range}.")
        logging.info("Config:power => ok !")

        return self.power_range
