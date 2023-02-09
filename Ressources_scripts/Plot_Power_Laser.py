import os
from pathlib import Path
import json
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import numpy as np

CUR_DIR = os.path.dirname(__file__)
DATA = Path(f"{CUR_DIR}")
DATA.mkdir(exist_ok=True)
DATA = DATA.parent
root = DATA/ "Ressources" / "dict_power_angle_conversion.json"



class plot:

    def __init__(self):

        pass

    def _read_power_angle_conversion(self, wavelength):
        with open(root, 'r') as f:
            dict = json.load(f)
        data = dict[f"{wavelength}"]
        return data

    def convert(self):
        x = []
        y = []
        z = []
        for wavelength in range(680, 1075, 5):
            data=self._read_power_angle_conversion(wavelength)
            for key, value in data.items():
                if int(key) > 30 and int(key) <125:
                    x.append(int(wavelength))
                    y.append(int(key))
                    z.append(float(1000*value))

        return np.array(x), np.array(y), np.array(z)

    def _get_data_from_angle(self, angle):
        x = []
        y = []
        for wavelength in range(680, 1081, 1):
            data=self._read_power_angle_conversion(wavelength)
            x.append(int(wavelength))
            intensity = data[f"{angle}"]
            y.append(int(intensity))
        return np.array(x), np.array(y)


    def graph(self, x, y, name):
        fig, ax = plt.subplots()
        plt.xlabel('Wavelength / nm')
        plt.ylabel('Laser intensity / cps')
        ax.plot(x, y, color='black')
        fig.savefig(DATA / f"{name}.png")

    def plot_3d(self, x, y, z, name):
        print("x = ", x)
        print("#"*100)
        print("y = ", y)
        print("#" * 100)
        print("z = ", z)
        print("#" * 100)
        fig, ax = plt.subplots()
        ax = plt.axes(projection='3d')
        ax.set_xlabel("Wavelength / nm")
        ax.set_ylabel("Angle / °")
        ax.set_zlabel("Laser intensity / cps")
        #ax.scatter(x, y, z, c=z, cmap='BrBG', linewidth=1)
        ax.plot_trisurf(x, y, z, cmap='viridis')
        #ax.plot_surface(x, y, z, cmap='jet')
        fig.savefig(DATA / f"{name}.png")
        #plt.draw()
        for angle in range(-150, 0):
            ax.view_init(elev=10, azim=angle, vertical_axis="z")
            plt.draw()
            plt.pause(.01)
        plt.show()

if __name__ == "__main__":
    m = plot()
    #x1, y1 = m._get_data_from_angle(80)
    #m.graph(x1, y1, "80")
    x, y, z = m.convert()
    m.plot_3d(x, y, z, "e")
