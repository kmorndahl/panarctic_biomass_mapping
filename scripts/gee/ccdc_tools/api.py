"""
@license
Copyright 2019 Boston University

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

   https://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

@fileoverview Main entry point for the CCDC API.
@author Eric L. Bullock, PhD

Python port of the original JavaScript api file from gee-ccdc-tools
(users/parevalo_bu/gee-ccdc-tools:ccdcUtilities/api).

Usage:
    import sys, os
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'parevalo_bu_gee_ccdc_tools', 'python'))
    from api import Classification, CCDC, Inputs, Dates, Change, Misc
"""

from . import classification as Classification
from . import ccdc as CCDC
from . import inputs as Inputs
from . import dates as Dates
from . import change as Change
from . import misc as Misc

__all__ = [
    'Classification',
    'CCDC',
    'Inputs',
    'Dates',
    'Change',
    'Misc',
]
