"""
Created on 5/1/2023
@author: Jonathan DANIEL, jonathan.daniel@u-bordeaux.fr
"""
from seabreeze.spectrometers import Spectrometer, list_devices
import numpy as np
import logging

logging.basicConfig(level=logging.INFO,
                    filename="PyTPEF.log",
                    filemode="a",
                    format='%(asctime)s - %(levelname)s - %(message)s')

class Maya():

    def __init__(self, integration_time):
        """
        Connects to the device
        :param integration_time: int
        """
        self.devices = list_devices()

        logging.info(f"Talk_Spectro:List of devices connected: {self.devices}.")

        try:
            self.spectro=Spectrometer.from_first_available()

            logging.info(f"Talk_Spectro:Selected device: {self.spectro}.")

            self.spectro.integration_time_micros(integration_time * 1000)

            logging.info(f"Talk_Spectro:Integration time: {integration_time}.")

        except:

            logging.critical("Talk_Spectro:No device found !")

    def get_spectrum(self):
        """
        This function reads data from the spectro
        :return: array, array
        """
        wavelength=self.spectro.wavelengths()
        intensities=self.spectro.intensities(correct_dark_counts=True, correct_nonlinearity=True)

        logging.info("Talk_Spectro:Data collected from spectro => ok !")

        return wavelength, intensities

    def _data_read_scan(self, number_of_scan, intensity_dark):
        """
        This function reads n times (number of scans) the spectro and correct it with dark noise.
        :param number_of_scan: int
        :param intensity_dark: int or array
        :return: array, array, array
        """
        if type(intensity_dark) == int:
            intensity_dark = np.zeros(shape=2068, dtype=int)
        for n in range(int(number_of_scan)):
            wavelength, intensities= self.get_spectrum()
            data = np.zeros(shape=2068, dtype=int)
            intensity = np.zeros(shape=2068, dtype=int)

            logging.info(f"Talk_Spectro:Data_n: {data}.")

            for i in range(len(intensities)):
                data[i]=data[i]+intensities[i]-intensity_dark[i]
                intensity[i]=intensity[i] + intensities[i]
                if type(intensity_dark) == int :
                    pass
                else :
                    if data[i] == 0:

                        logging.info(f"Talk_Spectro:Change {data[i]} to 1.")

                        data[i] =1


            logging.info(f"Talk_Spectro:type(Data): {type(data)}.")

        #logging(f'Talk_Spectro:type(Dataread): {type(data)}.')

        return wavelength, data, intensity

