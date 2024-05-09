"""
Stub module for downloading the latest analysis.py from https://github.com/lwa-project/stress_tests
"""

import os
import time
from urllib import request as urlrequest


# URL to download
MODULE_URL = 'https://raw.githubusercontent.com/lwa-project/stress_tests/main/analysis.py'

# Where it lives
MODULE_PATH = os.path.dirname(os.path.abspath(__file__))
MODULE_NAME = os.path.join(MODULE_PATH, '_analysis.py')
MODULE_ETAG = os.path.join(MODULE_PATH, '_analysis.etag')

# Maximum file age to accept without checking for a newer version
MAX_AGE_SEC = 86400


# Get the current file age and entity tag (if it exists)
age = 1e6
etag = ''
if os.path.exists(MODULE_NAME):
    age = time.time() - os.path.getmtime(MODULE_NAME)
    if os.path.exists(MODULE_ETAG):
        with open(MODULE_ETAG, 'r') as fh:
            etag = fh.read()
            
# If the file is more than MAX_AGE_SEC old, check for an update
if age > MAX_AGE_SEC:
    request = urlrequest.Request(MODULE_URL)
    opener = urlrequest.build_opener()
    data = opener.open(request)
    if data.headers['etag'] != etag:
        with open(MODULE_NAME, 'wb') as fh:
            fh.write(data.read())
        with open(MODULE_ETAG, 'w') as fh:
            fh.write(data.headers['etag'])
            
# Load in everything from the module
from _analysis import *
from _analysis import __all__
