import random
import numpy as np
import logging

logging.basicConfig(level=logging.INFO,
                    filename="PyTPEF.log",
                    filemode="a",
                    format='%(asctime)s - %(levelname)s - %(message)s')


class Maya():

    def __init__(self, integration_time):

        logging.info(f"Spectro_sim:integration_time: {integration_time}.")

    def get_spectrum(self):
        intensities=[]
        wavelength =[]
        for i in range (2000):
            wavelength.append(300.0+(i*0.4))
            i=random.randrange(0,i+1)
            a=random.randrange(1,10)
            value = (np.sin((i*(np.pi/180))/100)+np.sin((i*(np.pi/90))/100))/a
            intensities.append(value)
        wavelength=np.array(wavelength)
        intensities=np.array(intensities)
        return wavelength, intensities

    def _data_read_scan(self, number_of_scan, intensity_dark):

        logging.info(f"Spectro_sim:Number of scan= {number_of_scan}.")

        if type(intensity_dark) == int:
            intensity_dark = np.zeros(shape=2000, dtype=int)
        for n in range(int(number_of_scan)):
            wavelength, intensities= self.get_spectrum()
            data = np.zeros(shape=2000, dtype=int)
            intensity = np.zeros(shape=2000, dtype=int)
            logging.info(f"data: {data}.")
            for i in range(len(intensities)):
                data[i]=data[i]+intensities[i]-intensity_dark[i]
                intensity[i]=intensity[i] + intensities[i]
                if type(intensity_dark) == int :
                    pass
                else :
                    if data[i] == 0:
                        data[i] =1

            logging.info(f"Spectro_sim:type(data): {type(data)}.")
        logging.info(f'Spectro_sim:type(dataread) {type(data)}.')

        return wavelength, data, intensity
