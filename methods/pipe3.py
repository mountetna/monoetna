import numpy as np
from helper import helper

@helper
def func(a, b):
    return (np.array(a)*np.array(b)).tolist()