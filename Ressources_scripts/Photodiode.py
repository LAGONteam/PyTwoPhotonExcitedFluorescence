import Ressources_scripts.Control_PD as ctrl

#import Control_PD as ctrl
import numpy as np

with ctrl.CustomTLPM(0) as tlpm:
    pass

def open():
    with ctrl.CustomTLPM(0) as tlpm:
        a = tlpm.get_power()

    return True

def state():
    return True

def set_wavelength(wavelength):
    """
    Change the wavelength for photodiode correction
    """
    tlpm.wavelength = wavelength
    return tlpm.wavelength

def read():
    """
    Returns the data measured by the power meter
    """
    a = tlpm.get_power()
    if str(a) == "nan" or str(a) == "inf":
        return 0
    else:
        return a


def close():
    tlpm.close()
    return True

if __name__ == '__main__':

    for n in range(20):
        print("#"*1000)
        read()
        print("#"*1000)
    close()

