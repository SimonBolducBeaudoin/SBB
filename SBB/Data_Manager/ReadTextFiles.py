#!/bin/env/python
#! -*- coding: utf-8 -*-

"""Regroups function used to read text files"""

import numpy as _np
import glob as _glob
import re as _re

def extract_numbers_from_files(file_pattern, pattern=r'Computing : (\d+\.\d+) \[s\]'):
    numbers = []
    for file_path in _glob.glob(file_pattern):
        with open(file_path, 'r') as file:
            for line in file:
                match = re.search(pattern, line)
                if match:
                    numbers.append(float(match.group(1)))
    return _np.r_[numbers]