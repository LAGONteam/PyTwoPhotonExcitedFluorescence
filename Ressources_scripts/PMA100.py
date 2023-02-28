"""
Created on 5/1/2023
@author: Jonathan DANIEL, jonathan.daniel@u-bordeaux.fr
"""
from pyThorlabsPM100x.driver import ThorlabsPM100x

import logging

logging.basicConfig(level=logging.INFO,
                    filename="PyTPEF.log",
                    filemode="a",
                    format='%(asctime)s - %(levelname)s - %(message)s')

class PMA100():

    def __init__(self):
        super().__init__()
        self.powermeter = ThorlabsPM100x()
        #available_wavelength=(185, 500, 860, 900, 980, 1050, 1100, 1310, 1550)#Only for S310C thermal detector
        self.available_devices = self.powermeter.list_devices() #Check which devices are available
        print("Devices detected : ", self.available_devices)
        self.i = 0
        self.x= []
        self.y = []

        logging.info("PMA100:Initialization => ok !")

    def Connect(self):
        """
        This function connects the first available device.
        :return: bool
        """
        try :
            self.powermeter.connect_device(device_addr = self.available_devices[0][0]) #Connect to the first available device
            print(f"Successful connection to {self.available_devices[0][0]}")
            (maxWL, minWL) = self.powermeter.read_min_max_wavelength()  # read max and min available wavelengths
            print(f"The sensing range of the connected sensor is [{minWL} nm; {maxWL} nm].")
        except:

            logging.error("PMA100:No device found !")

            return False
        return True

    def Set_Power_Range(self, new_power_range):
        """
        This functions set the power range i.e. 1W, 100 mW, 10 mW, 1mW etc...
        :param new_power_range: int
        :return: state of a function
        """

        logging.info("PMA100:Starting Set_Power_Range.")
        logging.info(f"PMA100:new_power_range: {new_power_range}.")

        if new_power_range > 0:
            self.powermeter.move_to_next_power_range(direction=+1)
        elif new_power_range < 0:
            self.powermeter.move_to_next_power_range(direction=-1)

        logging.info(f"PMA100:The current power range is: {self.powermeter.power_range}.")  # print current power range
        logging.info("PMA100:Set_Power_Range => ok !")

        return self.powermeter.power_range

    def Read_Power_Range(self):
        """
        This function reads the power range.
        :return: state of a function
        """

        logging.info(f"PMA100:The current power range is: {self.powermeter.power_range}.")

        return self.powermeter.power_range

    def Set_Wavelength(self, new_wavelength):
        """
        This function set the wavelength
        :param new_wavelength: int
        :return: None
        """
        self.powermeter.wavelength = new_wavelength

        logging.info(f"PMA100:The powermeter is now set to {self.powermeter.wavelength} nm.")

    def Set_Zero(self):
        """
        This function set the minimal value of power to 0.
        :return: None
        """
        self.powermeter.set_zero()

        logging.info("PMA100:Zero is done.")

    def Read_Power(self):
        """
        This function reads the current power value.
        :return: lst, lst
        """
        self.i+=1

        logging.info(f"PMA100:x= {self.i}.")

        self.x=self.i
        self.y= (self.powermeter.power[0])*10*1000 #because the bs is 80/20, *1000 is to set W to mW
        self.y=round(self.y, 1)
        logging.info(f"PMA100:y= {self.y}.")

        return self.x, self.y

    def Clear(self):
        """
        This function clear the lists x and y.
        :return: None
        """

        logging.info("PMA100:Clear data.")

        self.x.clear()
        self.y.clear()

    def Disconnect(self):
        """
        This function disconnects the device.
        :return: bool
        """
        try :
            self.powermeter.disconnect_device() #disconnect the device

            logging.info(f"PMA100:Successful disconnection of {self.available_devices[0][0]}.")

        except:

            logging.warning(f"PMA100:Unable to disconnect {self.available_devices[0][0]}.")

            return False
        return True

    def Power_Units(self):
        """
        This function return the unit of the power.
        :return: state of a function
        """

        logging.info(f"PMA100:self.powermeter.power_units: {self.powermeter.power_units}.")

        return self.powermeter.power_units

