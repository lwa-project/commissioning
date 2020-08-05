"""
Stub module for getting the latest version of DRX/HDF5/data.py
"""

from __future__ import print_function

import os
import shutil

# Make sure we have the newest version of the file.  If not, copy it over
copy = False
if os.path.exists('_data.py'):
    if os.path.getmtime('../DRX/HDF5/data.py') > os.path.getmtime('_data.py'):
        copy = True
else:
    copy = True
    
if copy:
    shutil.copy('../DRX/HDF5/data.py', '_data.py')
    
# Load in everything from the module
from _data import *

