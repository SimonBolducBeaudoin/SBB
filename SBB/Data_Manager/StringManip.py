#!/bin/env/python
#! -*- coding: utf-8 -*-

"""Regroups function used to modifie arrays of strings"""

import numpy as _np

def remove_trailing(arr, separator=['-', '&'], trail='0'):
    """
    Remove the last occurrence of the specified trailing character from each number in the given array of strings.
    
    # Example of usage:
    arr = np.array(['100-200&3000', '400-500&6000', '700-800&9000'])
    modified_arr = remove_one_trailing_zero(arr, separator=['-', '&'], trail='0')
    >>> print(modified_arr)
    ['10-20&300' '40-50&600' '70-80&900']


    Parameters:
    arr (numpy.ndarray): An array of strings containing numbers separated by specified separators.
    separator (list, optional): A list containing separator characters used to split the strings.
                                Defaults to ['-','&'].
    trail (str, optional): The trailing character to remove from each number.
                           Defaults to '0'.

    Returns:
    numpy.ndarray: An array of strings where the last occurrence of the specified trailing character in each number is removed.
    """
    modified_arr = _np.empty_like(arr, dtype=arr.dtype)

    for i, label in enumerate(arr):
        parts = label.split(separator[0])
        modified_parts = []
        for part in parts:
            numbers = part.split(separator[1])
            modified_numbers = []
            for number in numbers:
                if trail:
                    modified_number = number[::-1].replace(trail, '', 1)[::-1]
                else:
                    modified_number = number.rstrip(trail)
                modified_numbers.append(modified_number)
            modified_parts.append(separator[1].join(modified_numbers))
        modified_label = separator[0].join(modified_parts)
        modified_arr[i] = modified_label
    return modified_arr