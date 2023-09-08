"""
Stub module for getting the latest version of DRX/HDF5/data.py
"""

import os
import shutil

MODULE_PATH = os.path.dirname(os.path.abspath(__file__))
MODULE_NAME = os.path.join(MODULE_PATH, '_data.py')
SOURCE_NAME = os.path.join(MODULE_PATH, '../DRX/HDF5/data.py')

# Make sure we have the newest version of the file.  If not, copy it over
copy = False
if os.path.exists(MODULE_NAME):
    if os.path.getmtime(SOURCE_NAME) > os.path.getmtime(MODULE_NAME):
        copy = True
else:
    copy = True
    
if copy:
    shutil.copy(SOURCE_NAME, MODULE_NAME)
    
# Load in everything from the module
from _data import *
