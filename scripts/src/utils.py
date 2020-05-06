#!/usr/bin/env python3
"""File containing versatile utility funcitions that may be used
for different scripts or functions.
"""

def dict_diff(a, b):
    """Returns a dictionary result of the elements only present in one
    of the inputs
    """
    return {k: a[k] if k in a else b[k] for k in a.keys() ^ b.keys()}


## END