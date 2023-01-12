import logging

logging.basicConfig(level=logging.INFO,
                    filename="PyTPEF.log",
                    filemode="a",
                    format='%(asctime)s - %(levelname)s - %(message)s')

class Task():

    def __init__(self):
        pass

    @property
    def ai_channels(self, tr):
        """
        :class:`nidaqmx._task_modules.ai_channel_collection.AIChannelCollection`:
            Gets the collection of analog input channels for this task.
        """
        logging.info(f"Photodiode_sim:Task.ai_channels tr: {tr}.")

        add_ai_voltage_chan(tr)
        return


class add_ai_voltage_chan():

    def __init__(self, value):
        self.value =value

        logging.info(f"Photodiode_sim:Task.ai_channels self.value: {self.value}.")