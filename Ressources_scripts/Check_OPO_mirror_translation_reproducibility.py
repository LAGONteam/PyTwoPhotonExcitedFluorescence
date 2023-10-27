import talk_to_elliptec_devices as motor
print("MOTORS OK")
import talk_to_spectro as spectro
print("SPECTRO OK")
import coherent_laser as laser
print("LASER OK")
import numpy as np

from pathlib import Path
import os
import matplotlib.pyplot as plt

INTEGRATION_TIME = 100
ANGLE_VALUES = [0, 5, 10, 20, 30, 35, 40, 45, 60, 90]
WAVELENGTH_RANGE = [680, 730, 750, 775, 800, 900, 1000, 1080]
I = 0
THETA = np.array(ANGLE_VALUES)
CUR_DIR = os.path.dirname(__file__)
ROOT_FOLDER = Path(f"{CUR_DIR}/Metrology")
ROOT_FOLDER.mkdir(exist_ok=True)

class Measure:

    def __init__(self):
        print("Start")
        laser.Chameleon().openShutterBlocking()
        self.move_opo_mirror(False)
        self.full_area = np.zeros((8,10))
        self.full_area_opo = np.zeros((8, 10))
        self.full_area_back_from_position = np.zeros((8,10))

    def busy(self, wavelengths,back_from_position, OPO):
        i = 0
        for angle in ANGLE_VALUES:
            print(f"{OPO}")
            print(angle)
            wavelength, fluorescence = self.measure(angle)
            power = self.measure_power()
            self.aggregate_data(wavelength, fluorescence, power, i)
            i += 1
        self.calculate_area(back_from_position, OPO)
        self.save_data(wavelengths, OPO)
        return

    def run(self):
        count =0
        for wavelengths in WAVELENGTH_RANGE:
            self.clean_table()
            motor.Power_detector().change_wavelength(int(wavelengths))
            print(wavelengths)
            laser.Chameleon().setWavelengthBlocking(wavelengths)
            self.control_shutter_laser(True)
            self.busy(wavelengths, False, False)
            self.full_area[count]=np.array(self.area)
            self.clean_table()
            self.move_opo_mirror(True)
            self.busy(wavelengths, False, True)
            self.full_area_opo[count]=np.array(self.area_opo)
            self.clean_table()
            self.move_opo_mirror(False)
            self.busy(wavelengths, True, False)
            self.control_shutter_laser(False)
            self.full_area_back_from_position[count]=np.array(self.area_back)
            count+=1
        self.plot()
        self.finish()

    def move_opo_mirror(self, opo_position):
        self.control_shutter_laser(False)
        motor.Linear_stage().pumping_opo(opo_position)
        self.control_shutter_laser(True)
        return

    def tune_laser_power(self, angle):
        motor.RotationMount().spin_to_position(angle)
        return

    def measure_fluorescence(self):
        wavelength, fluorescence = spectro.Maya(INTEGRATION_TIME).get_spectrum()
        return wavelength, fluorescence

    def measure_power(self):
        self.power = []
        for n in range(100):
            self.power.append(1000*(motor.Power_detector().read()))
            # self.power.append(n)
            # sleep(0.1)
        return np.array(self.power)

    def control_shutter_laser(self, open):
        if open:
            print("open")
            motor.Optical_Shutter().open()
        else:
            print("close")
            motor.Optical_Shutter().close()
        return

    def measure(self, angle):
        self.tune_laser_power(angle)
        wavelength, fluorescence = self.measure_fluorescence()
        return wavelength, fluorescence

    def clean_table(self):
        self.fluorescence_table = np.zeros((10, 2068))
        self.power_table = np.zeros((10, 100))

    def aggregate_data(self, wavelength, fluorescence, power, i):
        global I
        if I ==0:
            self.wavelength_table = wavelength

            I = 1
        else:
            self.fluorescence_table[i] = np.array(fluorescence)
            self.power_table[i] = np.array(power)
        return

    def calculate_area(self, back_from_position, OPO):
        i = 0
        area = []
        for n in self.fluorescence_table:
            temp = np.trapz(self.fluorescence_table[i], self.wavelength_table)
            area.append(temp)
            i+=1
        if OPO:
            self.area_opo = np.array(area)
        else:
            if back_from_position:
                self.area_back = np.array(area)
            else:
                self.area = np.array(area)
        return

    def plot(self):
        count = 0
        for n in WAVELENGTH_RANGE:
            fig, ax = plt.subplots()
            ax.plot(THETA, self.full_area[count], color='black')
            ax.plot(THETA, self.full_area_opo[count], color='red')
            ax.plot(THETA, self.full_area_back_from_position[count], color='blue')
            count+=1
            plt.xlabel('Angle / °')
            plt.ylabel('Fluorescence intensity / cps')
            fig.savefig(ROOT_FOLDER / f"{n}.png")

    def save_data(self, wavelengths, OPO):
        root_file = ROOT_FOLDER / f"OPO_translation_reproducibility_pumping_OPO_{OPO}_{wavelengths}_fluorescence.txt"
        root_file.touch()
        i = 0
        head_data = f"Wavelength (nm) \t Angle 0 (°) \t Angle 5 (°) \t Angle 10 (°) \t Angle 20 (°) \t Angle 30 (°) \t Angle 35 (°) \t Angle 40 (°) \t Angle 45 (°) \t Angle 60 (°) \t Angle 90 (°) \n"
        with open(root_file, "w") as w:
            w.write(head_data)
        for n in self.wavelength_table:
            data = f"{self.wavelength_table[i]} \t {self.fluorescence_table[0][i]} \t {self.fluorescence_table[1][i]} \t {self.fluorescence_table[2][i]} \t {self.fluorescence_table[3][i]} \t {self.fluorescence_table[4][i]} \t {self.fluorescence_table[5][i]} \t {self.fluorescence_table[6][i]} \t {self.fluorescence_table[7][i]} \t {self.fluorescence_table[8][i]} \t {self.fluorescence_table[9][i]} \n"
            with open(root_file, "a") as a:
                a.write(data)
            i+=1
        root_file = ROOT_FOLDER / f"OPO_translation_reproducibility_{wavelengths}_power.txt"
        root_file.touch()
        head_data = f"Angle 0 (mW) \t Angle 5 (mW) \t Angle 10 (mW) \t Angle 20 (mW) \t Angle 30 (mW) \t Angle 35 (mW) \t Angle 40 (mW) \t Angle 45 (mW) \t Angle 60 (mW) \t Angle 90 (mW) \n"
        i=0
        with open(root_file, "w") as w:
            w.write(head_data)
        for n in self.power_table[0]:
            data = f"{self.power_table[0][i]} \t {self.power_table[1][i]} \t {self.power_table[2][i]} \t {self.power_table[3][i]} \t {self.power_table[4][i]} \t {self.power_table[5][i]} \t {self.power_table[6][i]} \t {self.power_table[7][i]} \t {self.power_table[8][i]} \t {self.power_table[9][i]} \n"
            with open(root_file, "a") as a:
                a.write(data)
            i+=1
        return

    def finish(self):
        laser.Chameleon().closeShutterBlocking()
        motor.Power_detector().close()
        spectro.Maya(INTEGRATION_TIME).shut_down()
        return

if __name__ == '__main__':

    a = Measure()
    a.run()
