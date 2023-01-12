"""
Created on Mon Oct  7 23:44:38 2019
@author: williamstoy
"""

import serial
import time
import re
import logging

logging.basicConfig(level=logging.INFO,
                    filename="PyTPEF.log",
                    filemode="a",
                    format='%(asctime)s - %(levelname)s - %(message)s')

class Chameleon:
    ##########################################################################
    ## INITIALIZATION FUNCTIONS
    ##########################################################################
    def __init__(self, com_port='COM5', verbose=False):
        try :
            self.verbose = verbose
            self.ser = serial.Serial(com_port, 19200, timeout=1)  # open serial port
            # initialize the connection by writing a carriage return
            self.sendCommand(b'')
            time.sleep(0.05)
            self.current_wavelength = 800

            logging.info("Laser:Laser initialized => ok !")

        except:

            logging.error("Laser:Laser cannot be initialized !")


    ##########################################################################
    ## COMMAND / RESPONSE PRIMITIVE FUNCTIONS
    ##########################################################################
    def sendCommand(self, cmd):

        logging.info(f"Laser:sendCommand:cmd: {cmd}.")

        self.ser.write(cmd + b'\r')

    def getResponse(self, cmd):

        logging.info(f"Laser:getResponse:cmd: {cmd}.")

        response = self.ser.read(50)
        # format the command so that it can be taken out with a regular expression
        cmd = re.sub(b'[\?]', b'\?', cmd)
        cmd = re.sub(b'[\=]', b'\=', cmd)
        # strip out \r and \n, CHAMELEON> and the command
        response = re.sub(b'[\\r|\\n|(CHAMELEON\>)|(' + cmd + b')|\s]', b'', response)

        if (self.verbose):

            logging.info(f"Laser:response: {response}.")

        return response

    def sendCmdGetResponse(self, cmd):

        logging.info(f"Laser:sendCmdGetResponse:cmd: {cmd}.")

        self.sendCommand(cmd)
        time.sleep(0.05)
        return self.getResponse(cmd)

    ##########################################################################
    ## WAVELENGTH FUNCTIONS
    ##########################################################################
    def setWavelength(self, n):

        logging.info(f"Laser:setWavelength:n: {n}.")

        self.current_wavelength = n
        if (self.verbose):

            logging.info(f'Laser:Tuning the laser to {n} nm')

        self.sendCmdGetResponse(b'vw=' + bytes(str(n), 'utf-8'))

    def queryTunedStatus(self):
        return self.sendCmdGetResponse(b'?ts')

    def setWavelengthBlocking(self, n):

        logging.info(f"Laser:setWavelengthBlocking:n: {n}.")

        self.setWavelength(n)
        while (self.queryTunedStatus() != b'0'):
            if (self.verbose):

                logging.info('Laser:Waiting for tuning to finish...')

            time.sleep(0.25)
        if (self.verbose):

            logging.info('Laser:Laser is Tuned')

    ##########################################################################
    ## SHUTTER FUNCTIONS
    ##########################################################################
    def openShutter(self):
        if (self.verbose):

            logging.info('Laser:Opening the Shutter')

        self.sendCmdGetResponse(b's=1')

    def closeShutter(self):
        if (self.verbose):

            logging.info('Laser:Closing the Shutter')

        self.sendCmdGetResponse(b's=0')

    def queryShutterStatus(self):
        return self.sendCmdGetResponse(b'?s')

    def openShutterBlocking(self):
        self.openShutter()
        while (self.queryShutterStatus() != b'1'):
            if (self.verbose):

                logging.info('Laser:Waiting for the shutter to open..')

            time.sleep(0.25)
        if (self.verbose):

            logging.info('Laser:Shutter is OPEN')

    def closeShutterBlocking(self):
        self.closeShutter()
        while (self.queryShutterStatus() != b'0'):
            if (self.verbose):

                logging.info('Laser:Waiting for the shutter to close..')

            time.sleep(0.25)
        if (self.verbose):

            logging.info('Laser:Shutter is CLOSED')

    def queryRelativeHumidity(self):
        if (self.verbose):

            logging.info('Laser:Getting Relative Humidity')

        self.sendCmdGetResponse(b'?rh')

    ##########################################################################
    ## SERIAL FUNCTIONS
    ##########################################################################
    def close(self):
        self.ser.close()
