"""
Module that is called at the begin of the simulation to set the log file trough the logging library
"""

__author__ = "Enrico Prataviera"
__credits__ = ["Enrico Prataviera"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Enrico Prataviera"


import os
import logging

from eureca_building.config import CONFIG

# Logging file
root_logger = logging.getLogger()
root_logger.setLevel(logging.ERROR)  # or whatever
handler = logging.FileHandler(
    os.path.join(CONFIG.output_path, "logging.log"), "w", "utf-8"
)  # or whatever
formatter = logging.Formatter(
    "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)  # or whatever
handler.setFormatter(formatter)  # Pass handler as a parameter, not assign
root_logger.addHandler(handler)
