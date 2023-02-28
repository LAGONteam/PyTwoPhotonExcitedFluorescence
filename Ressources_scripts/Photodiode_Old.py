from ThorlabsPM100 import ThorlabsPM100
import pyvisa
import time

rm = pyvisa.ResourceManager()
inst = rm.open_resource('USB0::0x1313::0x807A::M00795996::INSTR')
POWER_METER = ThorlabsPM100(inst=inst)
POWER_METER.sense.correction.collect.zero.initiate
print("Power range limit:", POWER_METER.sense.power.dc.range.upper)
print("Power range limit:", POWER_METER.sense.power.dc.range.maximum_upper)
POWER_METER.sense.power.dc.range.auto='ON'

def state():
    """
    Check if the power meter is connected
    """
    return True

def read():
    """
    Returns the data measured by the power meter
    """
    return POWER_METER.read

def zero():
    """
    Sets the dark offset
    """
    try:
        POWER_METER.sense.correction.collect.zero.initiate
        return True
    except:
        return False

def get_wavelength():
    """
    Change the wavelength for photodiode correction
    """
    return POWER_METER.sense.correction.wavelength

def set_wavelength(wavelength):
    """
    Change the wavelength for photodiode correction
    """
    POWER_METER.sense.correction.wavelength = wavelength
    return POWER_METER.sense.correction.wavelength

def get_configure():
    """
    Returns the actual mode of the power meter
    """
    return POWER_METER.getconfigure

if __name__ == '__main__':
    print(read())
    print(a)
    print(get_wavelength())
    print(set_wavelength(690))
    time.sleep(2)
    a = get_configure()