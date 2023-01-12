import serial
import time
import re
import logging

import logging

logging.basicConfig(level=logging.INFO,
                    filename="PyTPEF.log",
                    filemode="a",
                    format='%(asctime)s - %(levelname)s - %(message)s')

class Chameleon:
    ##########################################################################
    ## INITIALIZATION FUNCTIONS
    ##########################################################################
    def __init__(self, verbose=False):
        self.verbose = verbose
        self.current_wavelength = 800

    ##########################################################################
    ## COMMAND / RESPONSE PRIMITIVE FUNCTIONS
    ##########################################################################
    def sendCommand(self, cmd):

       logging.info(f"Laser_sim:LASER:cmd: {cmd}.")


    def getResponse(self, cmd):

        if (self.verbose):
            response = True

            logging.info(f"Laser_sim:LASER:response: {response}.")
            logging.info(f"Laser_sim:LASER:cmd: {cmd}.")

        return response

    def sendCmdGetResponse(self, cmd):
        self.sendCommand(cmd)
        time.sleep(0.05)
        return self.getResponse(cmd)

    ##########################################################################
    ## WAVELENGTH FUNCTIONS
    ##########################################################################
    def setWavelength(self, n):
        self.current_wavelength = n
        return True

    def queryTunedStatus(self):
        return True

    def setWavelengthBlocking(self, n):
        self.setWavelength(n)

        logging.info(f"Laser_sim:LASER:Wavelength tuning to {n}")

    ##########################################################################
    ## SHUTTER FUNCTIONS
    ##########################################################################
    def openShutter(self):

        logging.info('Laser_sim:LASER:Opening the Shutter')


    def closeShutter(self):

        logging.info('Laser_sim:LASER:Closing the Shutter')


    def queryShutterStatus(self):
        return True

    def openShutterBlocking(self):
        self.openShutter()


    def closeShutterBlocking(self):
        self.closeShutter()
        
        logging.info('Laser_sim:LASER:Shutter is CLOSED')

    def queryRelativeHumidity(self):
        pass

    ##########################################################################
    ## SERIAL FUNCTIONS
    ##########################################################################
    def close(self):

        logging.info("Laser_sim:LASER:Close Laser")
