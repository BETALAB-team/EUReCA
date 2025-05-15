"""
Initializes the logging system for the simulation.

This module configures the root logger to write error messages to a file located
in CONFIG.output_path/logging.log using UTF-8 encoding and a standard formatter.
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
