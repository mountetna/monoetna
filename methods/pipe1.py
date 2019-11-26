import numpy as np
from helper import helper

@helper
def func(arr1, arr2):
    return (np.array(arr1) + np.array(arr2)).tolist()