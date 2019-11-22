import numpy as np

def func(arr):
    if len(np.shape(arr)) == 1:
        arr = np.reshape(arr, [len(arr),1])
        
    return np.dot(arr, arr.T).tolist(), np.dot(arr.T, arr).tolist()