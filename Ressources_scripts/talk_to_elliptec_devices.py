"""
Created on 5/1/2023
@author: Jonathan DANIEL, jonathan.daniel@u-bordeaux.fr
"""
import time
import numpy as np
import elliptec
#from Ressources_scripts import PMA100
import os
from pathlib import Path
import json
#from coherent_laser import Chameleon
from Ressources_scripts.coherent_laser import Chameleon
import logging
from Ressources_scripts import Photodiode
#import Photodiode

logging.basicConfig(level=logging.INFO,
                    filename="PyTPEF.log",
                    filemode="a",
                    format='%(asctime)s - %(levelname)s - %(message)s')

CUR_DIR = os.path.dirname(__file__)
DATA = Path(f"{CUR_DIR}/Datas")
DATA.mkdir(exist_ok=True)
print(DATA)
root = Path(f"{CUR_DIR}").parent
root = root / "Ressources" / "dict_power_angle_conversion.json"
CONTROLLER= elliptec.Controller('COM4')
ROTATION= elliptec.Rotator(CONTROLLER)
ROTATION.home()



# from ThorlabsPM100 import ThorlabsPM100
# import pyvisa
#
# rm = pyvisa.ResourceManager()
# inst = rm.open_resource('USB0::0x1313::0x807A::M00795996::INSTR')
# power_meter_USB=ThorlabsPM100(inst=inst)
# power_meter_USB.configure.scalar.power
# power_meter_USB.sense.correction.collect.zero.initiate



#POWER_METER= PMA100.PMA100()
#POWER_METER.Connect()


class power_angle_conversion:

    def __init__(self):
        self.dict_power_angle_conversion = {}

        logging.info("Talk_elliptec:Initialisation of power_angle_conversion => ok !")

    def power_angle_conversion(self, wavelength):
        """
        This function read the power of the PMA100 at a given half-wave plate angle and updates the dictionnary with
        the values at a given wavelength
        :param wavelength: int
        :return: None
        """

        logging.info("Talk_elliptec:Starting power_angle_conversion.")

        self.dict_power_angle_conversion.update({wavelength: {0: ""}})
        for angle in range(90):
            if angle%2 == 0:

                logging.info(f"Talk_elliptec:The angle is {angle}.")

                RotationMount().rotation.set_angle(angle)
                time.sleep(0.2) #for power stabilization

                logging.info(f"Talk_elliptec:The real angle is: {RotationMount().get_angle()}.")

                # power = power_meter_USB.read
                power  = Photodiode.read()

                logging.info(f"Talk_elliptec:The corresponding power is {power} mW.")

                self.dict_power_angle_conversion[wavelength].update({angle: power})

        logging.info("power_angle_conversion => ok !")

    def measure_power_angle_conversion(self):
        """
        This function measures and saves the power of the laser for all wavelengths and all half-wave plate angle.
        :return: None
        """

        logging.info("Talk_elliptec:Starting measure_power_angle_conversion.")

        for wavelength in range(680, 1081, 1):
            # power_meter_USB.sense.correction.wavelength=int(wavelength)
            Photodiode.set_wavelength(wavelength)

            logging.info(f"Talk_elliptec:The wavelength is {wavelength} nm.")

            Chameleon().setWavelengthBlocking(wavelength)

            logging.info(f"Talk_elliptec:Laser move to {wavelength} => ok !")

            Chameleon().openShutterBlocking()

            logging.info("Talk_elliptec:Laser shutter is open !")

            time.sleep(10)
            self.power_angle_conversion(wavelength)
            Chameleon().closeShutterBlocking()

            logging.info("Talk_elliptec:Laser shutter is close.")

        time.sleep(0.1)
        with open(root, "w") as g:
            json.dump(self.dict_power_angle_conversion, g)

        logging.info("Talk_elliptec:measure_power_angle_conversion => ok !")

    def read_power_angle_conversion(self, wavelength):
        """
        Unused function
        """
        """
        This function reads the power angle conversion dictionnary and returns the data corresponding to the wavelength
        :param wavelength:  int
        :return: dict
        """

        logging.info("Talk_elliptec:Starting read_power_angle_conversion.")
        logging.info(f"Talk_elliptec:wavelength: {wavelength}.")

        with open(root, 'r') as f:
             dict = json.load(f)
        data = dict[f"{wavelength}"]

        logging.info("Talk_elliptec:read_power_angle_conversion => ok !")

        return data

    def update_power_angle_conversion(self, wavelength, angle, power):

        logging.info("Talk_elliptec:Starting update_power_angle_conversion.")
        logging.info(f"Talk_elliptec:wavelength: {wavelength}.")
        logging.info(f"Talk_elliptec:angle: {angle}.")
        logging.info(f"Talk_elliptec:power: {power}.")


        with open(root, 'r') as f:
             dict = json.load(f)
        try:
            print(f"Check old entry dict[{wavelength}][{angle}]: ", dict[f"{wavelength}"][f"{angle}"])
        except KeyError:
            print("The angle is created")

        dict[f"{wavelength}"].update({f"{angle}": power})


        print(f"Check new entry dict[{wavelength}][{angle}]: ", dict[f"{wavelength}"][f"{angle}"])
        with open(root, "w") as w:
            json.dump(dict, w)



class RotationMount:

    """Initialization of the Rotation Mount"""
    """ 
    Important steps to do !
    #step 1 => get home offset
    #step2 => perform homing 
    """

    def __init__(self):
        self.controller= CONTROLLER
        self.rotation= ROTATION
        """self.controller = elliptec.Controller('COM4')
        self.rotation=elliptec.Rotator(self.controller)
        info = self.rotation.get('info')
        self.rotation.home()"""
        time.sleep(1)
        self.dico_power = {}
        self.min = False
        self.max = False
        self.dict_power_angle_conversion = {}
        #self.p = POWER_METER
        #self.p.Connect()

        logging.info("Talk_elliptec:Initialisation of RotationMount => ok !")

    def get_angle(self):
        """
        This function finds at which angle (in degrees) the rotator is at the moment.
        :return: int
        """

        logging.info("Talk_elliptec:Starting get_angle.")

        angle = elliptec.rotator.Rotator(controller=CONTROLLER).get_angle()

        logging.info(f"Talk_elliptec:angle: {angle}.")
        logging.info("Talk_elliptec:get_angle => ok !")

        return angle

    def read_power_meter(self):
        """
        This function reads the power of the PMA100 and means the values.
        :return: float
        """

        logging.info("Talk_elliptec:Starting read_power_meter.")

        x, y = self.p.Read_Power()
        table = np.mean(np.array(y))

        logging.info("Talk_elliptec:read_power_meter => ok !")

        return table

    def determine_angle_for_power(self, power_max, power_min, number_of_measure, wavelength):
        """
        This function determines to which angles, the half-wave plate must be sets in order to get desired laser power.
        :param power_max: int
        :param power_min: int
        :param number_of_measure: int
        :param wavelength: int
        :return: int, int, dict
        """

        logging.info("Talk_elliptec:Starting determine_angle_for_power.")
        logging.info(f"Talk_elliptec:power_max: {power_max}.")
        logging.info(f"Talk_elliptec:power_min: {power_min}.")
        logging.info(f"Talk_elliptec:number_of_measure: {number_of_measure}.")
        logging.info(f"Talk_elliptec:wavelength: {wavelength}.")

        self.power_max = float(power_max/1000)
        self.power_min = float(power_min/1000)
        self.number_of_measure = number_of_measure

        logging.info("Talk_elliptec:Start Angle for power ")

        '''Adapt the power range of the PM100A depending on the value of power to reach '''
        """if self.power_min > 400:
            p = self.p.Read_Power_Range()
            if p < 400:
                l =self.p.Set_Power_Range(+1)"""

        '''measure the power at each angle (0-180) to get the closest angle values for min and max power to reach'''
        """for angle in range(90):
            self.rotation.set_angle(angle)
            time.sleep(0.1)
            m =self.read_power_meter()*1000 #1000 is due to user input in mW while PM100A deal with W"""
        angle_power_measurement = self.read_power_angle_conversion(wavelength)
        for angle, power in angle_power_measurement.items():
            if power == self.power_min:
                self.angle_min = angle
                self.min = True
            elif power == self.power_max:
                self.angle_max = angle
                self.max = True
            else:
                self.dico_power[angle] = power
            if self.min == True and self.max == True :

                logging.info("Talk_elliptec:determine_angle_for_power stop at first level !")

                return
        if self.min == False:
            self.angle_min = self.min_angle()

            logging.info(f"Talk_elliptec:The minimum angle is: {self.angle_min}, for a power of : {self.dico_power[self.angle_min]} mW.")

            self.min = True
        if self.max == False:
            self.angle_max = self.max_angle()

            logging.info(f"Talk_elliptec:The minimum angle is: {self.angle_max}, for a power of : {self.dico_power[self.angle_max]} mW.")

            self.max = True
        if self.max == True and self.min == True:

            logging.info("Talk_elliptec:determine_angle_for_power => ok !")

            return self.angle_min, self.angle_max, self.dico_power

    def max_angle(self):
        """
        This function determines the closest angle to max power.
        :return: int
        """

        logging.info("Talk_elliptec:Starting max_angle.")

        max = self.closest_value(self.power_max, self.dico_power)
        for key, value in self.dico_power.items():
            if value == max :

                logging.info(f"Talk_elliptec:Key: {key}.")
                logging.info("Talk_elliptec:max_angle => ok !")

                return key

    def min_angle(self):
        """
        This function determines the closest angle to min power.
        :return: int
        """

        logging.info("Talk_elliptec:Starting min_angle.")

        min = self.closest_value(self.power_min, self.dico_power)
        for key, value in self.dico_power.items():
            if value == min :

                logging.info(f"Talk_elliptec:Key: {key}.")
                logging.info("Talk_elliptec:min_angle => ok !")

                return key

    def spin_to_position(self, position):
        """
        This function set the rotative motor (half-wave plate) to the angle input.
        :param position: int
        :return: None
        """

        logging.info("Talk_elliptec:Starting spin_to_position.")
        logging.info(f"Talk_elliptec:position: {position}.")

        self.rotation.set_angle(int(position))

        logging.info("Talk_elliptec:spin_to_position => ok !")

    def closest_value(self, value_to_reach, dictionnary):
        """
        This function returns the value of the input dict the closest the the value to rach input.
        :param value_to_reach: int
        :param dictionnary: dict
        :return: int
        """

        logging.info("Talk_elliptec:Starting closest_value.")
        logging.info(f"Talk_elliptec:value_to_reach: {value_to_reach}.")
        logging.info(f"Talk_elliptec:type(dictionnary): {type(dictionnary)}.")

        higher = []
        lower = []
        for key, value in dictionnary.items():
            if value > value_to_reach:
                #print(f"{value} > {value_to_reach}")
                higher.append(value)
            elif value < value_to_reach :
                #print(f"{value} < {value_to_reach}")
                lower.append(value)
        if len(higher) == 0:

            logging.info("Talk_elliptec:closest_value => ok !")

            return max(lower)
        if len(lower) == 0:

            logging.info("Talk_elliptec:closest_value => ok !")

            return min(higher)
        up = abs(min(higher)-value_to_reach)

        logging.info(f"Talk_elliptec:up: {up}.")

        down = abs(max(lower)-value_to_reach)

        logging.info(f"Talk_elliptec:down: {down}.")

        if up > down :

            logging.info("Talk_elliptec:closest_value => ok !")

            return max(lower)
        else :

            logging.info("Talk_elliptec:closest_value => ok !")

            return min(higher)

    def angles_for_measure(self, angle_min, angle_max, dictionary, number):
        """
        This function contructs the angle list to measure.
        :param angle_min: int
        :param angle_max: int
        :param dictionary: dict
        :param number: int
        :return: lst, lst, dict
        """

        logging.info("Talk_elliptec:Starting angles_for_measure.")
        logging.info(f"Talk_elliptec:angle_min: {angle_min}.")
        logging.info(f"Talk_elliptec:angle_max: {angle_max}.")
        logging.info(f"Talk_elliptec:type(dictionary): {type(dictionary)}.")
        logging.info(f"Talk_elliptec:number: {number}.")

        self.angle_min= angle_min
        self.angle_max = angle_max
        self.dico_power = dictionary
        self.number_of_measure = number
        list_of_angles = []
        list_of_power = []
        list_of_angles.append(self.angle_min)
        list_of_power.append(self.dico_power[self.angle_min])
        interval = (self.dico_power[self.angle_max] - self.dico_power[self.angle_min]) / (self.number_of_measure + 1)
        for n in range(self.number_of_measure):
            angle_ref = (n + 1)* interval + self.dico_power[self.angle_min]
            real_value = self.closest_value(angle_ref, self.dico_power)
            for key, value in self.dico_power.items():
                if value == real_value:
                    list_of_angles.append(key)
                    list_of_power.append(value)
        list_of_angles.append(self.angle_max)
        list_of_power.append(self.dico_power[self.angle_max])

        logging.info("Talk_elliptec:angles_for_measure => ok !")

        return list_of_angles, list_of_power, self.dico_power

    def power_angle_conversion(self, wavelength):
        """
        This function is used by an unused function ! Check this !
        :param wavelength:
        :return:
        """

        logging.info("Talk_elliptec:Starting power_angle_conversion.")
        logging.info(f"Talk_elliptec:wavelength: {wavelength}.")

        self.dict_power_angle_conversion.update({wavelength: {0: ""}})
        for angle in range(181):

            logging.info(f"Talk_elliptec:The angle is {angle}.")

            self.rotation.set_angle(angle)
            time.sleep(4)  # for power stabilization
            power = self.read_power_meter() * 1000
            self.p.Clear()

            logging.info(f"Talk_elliptec:The corresponding power is {power} mW.")

            self.dict_power_angle_conversion[wavelength].update({angle: power})

        logging.info("Talk_elliptec:power_angle_conversion => ok !")

    def measure_power_angle_conversion(self):
        """
        Unused function
        """

        logging.info("Talk_elliptec:Starting measure_power_angle_conversion.")

        self.p.Set_Zero()
        time.sleep(5)
        Chameleon().openShutterBlocking()

        logging.info("Talk_elliptec:Laser shutter is open !")

        for wavelength in range(680, 1081, 5):

            logging.info(f"Talk_elliptec:The wavelength is {wavelength} nm.")

            Chameleon().setWavelengthBlocking(wavelength)
            self.power_angle_conversion(wavelength)
        Chameleon().closeShutterBlocking()

        logging.info("Talk_elliptec:The laser shutter is close.")

        Chameleon().setWavelengthBlocking(800)

        logging.info("Talk_elliptec:The laser is tuned to 800 nm !")

        with open(root, "w") as g:
            json.dump(self.dict_power_angle_conversion, g)

        logging.info("Talk_elliptec:measure_power_angle_conversion => ok !")

    def read_power_angle_conversion(self, wavelength):
        """
        This function returns the angle-power data at a specific wavelength
        :param wavelength: int
        :return: dict
        """

        logging.info("Talk_elliptec:Starting read_power_angle_conversion.")
        logging.info(f"Talk_elliptec:wavelength: {wavelength}.")

        with open(root, 'r') as f:
            dict = json.load(f)
        data = dict[f"{int(wavelength)}"]

        logging.info("Talk_elliptec:read_power_angle_conversion => ok !")

        return data


if __name__ == '__main__':

    m = power_angle_conversion()
    m.measure_power_angle_conversion()


