import random
import time
import numpy as np
import os
from pathlib import Path
import json
from Ressources_scripts.Laser_Simulation import Chameleon
import logging

logging.basicConfig(level=logging.INFO,
                    filename="PyTPEF.log",
                    filemode="a",
                    format='%(asctime)s - %(levelname)s - %(message)s')

CUR_DIR = os.path.dirname(__file__)
temp = Path(f"{CUR_DIR}").parent
DATA = Path(f"{CUR_DIR}/Datas")
DATA.mkdir(exist_ok=True)
print(DATA)
root = Path(f"{temp}/Ressources/dict_power_angle_conversion.json")

class power_angle_conversion:

    def __init__(self):
        self.dict_power_angle_conversion = {}

    def power_angle_conversion(self, wavelength):
        self.dict_power_angle_conversion.update({wavelength: {0: ""}})
        for angle in range(90):
            if angle%2 == 0:

                logging.info(f"Motor_sim:The angle is {angle}.")

                time.sleep(0.2)

                logging.info("Motor_sim:The real angle is :", angle)

                power = RotationMount().read_power_meter() * 1000
                power = round(power, 2)

                logging.info(f"Motor_sim:The corresponding power is {power} mW.")

                self.dict_power_angle_conversion[wavelength].update({angle: power})

    def measure_power_angle_conversion(self):
        time.sleep(5)

        logging.info("Motor_sim:run")

        for wavelength in range(680, 950, 10):

            logging.info(f"Motor_sim:The wavelength is {wavelength} nm.")

            Chameleon().setWavelengthBlocking(wavelength)
            Chameleon().openShutterBlocking()
            #time.sleep(10)
            self.power_angle_conversion(wavelength)
            Chameleon().closeShutterBlocking()
        time.sleep(0.1)
        with open(root, "w") as g:
            json.dump(self.dict_power_angle_conversion, g)

    def read_power_angle_conversion(self, wavelength):
        with open(root, 'r') as f:
             dict = json.load(f)
        data = dict[f"{wavelength}"]
        return data

class RotationMount:

    def __init__(self):
        time.sleep(1)
        self.dico_power = {}
        self.min = False
        self.max = False
        self.dict_power_angle_conversion = {}
    def get_angle(self):
        ''' Finds at which angle (in degrees) the rotator is at the moment. '''
        angle = random.randint(0,90)
        return angle


    def read_power_meter(self):
        y = [random.randint(-10,100)*i for i in range(500)]
        table = np.mean(np.array(y))
        return table

    def determine_angle_for_power(self, power_max, power_min, number_of_measure, wavelength):
        self.power_max = float(power_max)
        self.power_min = float(power_min)
        self.number_of_measure = number_of_measure

        logging.info("Motor_sim:Start Angle for power ")

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
                return
        if self.min == False:
            self.angle_min = self.min_angle()

            logging.info(f"MOTOR:The minimum angle is: {self.angle_min}, for a power of: {self.dico_power[self.angle_min]} mW.")

            self.min = True
        if self.max == False:
            self.angle_max = self.max_angle()

            logging.info(f"MOTOR:The minimum angle is: {self.angle_max}, for a power of: {self.dico_power[self.angle_max]} mW.")

            self.max = True
        if self.max == True and self.min == True:
            return self.angle_min, self.angle_max, self.dico_power

    def max_angle(self):
        """determine the angle the closest to power_max value"""
        max = self.closest_value(self.power_max, self.dico_power)
        for key, value in self.dico_power.items():
            if value == max :
                return key

    def min_angle(self):
        """determine the angle the closest to power_min value"""
        min = self.closest_value(self.power_min, self.dico_power)
        for key, value in self.dico_power.items():
            if value == min :
                return key

    def spin_to_position(self, position):
        """Go to the angle given as an argument"""

        logging.info("Motor_sim:position")

    def closest_value(self, value_to_reach, dictionnary):
        """This function return the value of the dictionnary the closest of the value_to_reach"""
        higher = []
        lower = []
        for key, value in dictionnary.items():
            if value > value_to_reach:
                higher.append(value)
            elif value < value_to_reach :
                lower.append(value)
        if len(higher) == 0:
            return max(lower)
        if len(lower) == 0:
            return min(higher)
        up = abs(min(higher)-value_to_reach)

        logging.info(f"Motor_sim:up: {up}.")

        down = abs(max(lower)-value_to_reach)

        logging.info(f"Motor_sim:down: {down}.")

        if up > down :
            return max(lower)
        else :
            return min(higher)

    def angles_for_measure(self, angle_min, angle_max, dictionary, number):
        """The construction of angle list is (min, min+interval*n, max)"""
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
        return list_of_angles, list_of_power, self.dico_power

    def power_angle_conversion(self, wavelength):
        self.dict_power_angle_conversion.update({wavelength: {0: ""}})
        for angle in range(181):

            logging.info(f"Motor_sim:The angle is {angle}.")

            time.sleep(4)  # for power stabilization
            power = self.read_power_meter() * 1000

            logging.info(f"Motor_sim:The corresponding power is {power} mW.")

            self.dict_power_angle_conversion[wavelength].update({angle: power})

    def measure_power_angle_conversion(self):
        time.sleep(5)
        Chameleon().openShutterBlocking()
        for wavelength in range(680, 1081, 5):

            logging.info(f"Motor_sim:The wavelength is {wavelength} nm.")

            Chameleon().setWavelengthBlocking(wavelength)
            self.power_angle_conversion(wavelength)
        Chameleon().closeShutterBlocking()
        Chameleon().setWavelengthBlocking(800)
        with open(root, "w") as g:
            json.dump(self.dict_power_angle_conversion, g)

    def read_power_angle_conversion(self, wavelength):
        with open(root, 'r') as f:
            dict = json.load(f)
        data = dict[f"{int(wavelength)}"]
        return data


