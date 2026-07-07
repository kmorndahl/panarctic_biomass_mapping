"""
Python port of parevalo_bu/gee-ccdc-tools ccdcUtilities package.

Modules:
    api            -- Main entry point; re-exports all sub-modules as named attributes
    ccdc           -- Utility functions for filtering and modifying CCD output
    change         -- Utility functions for exploring and processing change information
    classification -- Utility functions for classifying CCDC results
    dates          -- Date conversion functions for the CCDC algorithm
    inputs         -- Utility functions for getting Landsat/Sentinel inputs for CCDC
    misc           -- Miscellaneous GEE asset utilities
    stats          -- Statistical distribution and test functions

Typical usage:
    import sys, os
    sys.path.insert(0, '/path/to/parevalo_bu_gee_ccdc_tools/python')
    from api import Dates, CCDC, Classification, Inputs, Change, Misc

Or import individual modules directly:
    from . import dates
    from . import ccdc
"""

from . import dates
from . import stats
from . import misc
from . import ccdc
from . import change
from . import inputs
from . import classification
from . import api

__all__ = [
    'dates',
    'stats',
    'misc',
    'ccdc',
    'change',
    'inputs',
    'classification',
    'api',
]
