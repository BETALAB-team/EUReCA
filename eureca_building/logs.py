"""
function for logs
"""

__author__ = "Enrico Prataviera"
__credits__ = ["Enrico Prataviera"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Enrico Prataviera"


import os
import logging

# Logging file
root_logger = logging.getLogger()
root_logger.setLevel(logging.WARNING)  # or whatever
handler = logging.FileHandler(
    os.path.join(".", "logging.log"), "w", "utf-8"
)  # or whatever
formatter = logging.Formatter(
    "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)  # or whatever
handler.setFormatter(formatter)  # Pass handler as a parameter, not assign
root_logger.addHandler(handler)
