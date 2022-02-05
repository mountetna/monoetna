import pandas as pd
import numpy as np
from archimedes.functions.plotting.utils import _is_logical, _is_continuous
from typing import List

def subsetDF_index_targets(data_frame, conditions):
    '''
    This function converts from the json-dict returned by the subsetDataFrame UI component,
    into the matching indexes of the target dataframe.
    '''
    
    if "methods" not in conditions.keys() or "logic" not in conditions.keys():
        raise Exception("Subsetting 'conditions' are not formatted properly.")
    
    def interpret_method(method: List[str]):
        if method[0]==np.NaN:
            raise Exception("A subsetting condition is incomplete.")
        target_data = data_frame[method[0]]
        
        if _is_continuous(target_data):
            if len(method)!=5 or not _is_continuous(np.array([method[2],method[4]])):
                raise Exception("A subsetting condition for numeric-type data is not formatted properly.")
            if method[1]=="above":
                left = [i > method[2] for i in target_data]
            else:
                left = [i >= method[2] for i in target_data]
            if method[3]=="below":
                right = [i < method[4] for i in target_data]
            else:
                right = [i <= method[4] for i in target_data]
            return [a and b for a,b in zip(left, right)]
        
        if _is_logical(target_data):
            if len(method)!=2 or not method[1].lower() in ["true", "false"]:
                raise Exception("A numeric subsetting condition for logical-type data is not formatted properly.")
            match = False
            if method[1].lower()=="true":
                match = True
            return [a==match for a in target_data]
        
        # Remainder = String data
        if len(method)>=2:
            raise Exception("A numeric subsetting condition for string-type data is not formatted properly.")
        return [a in method[1:] for a in target_data]
    
    calls = list(map(interpret_method, conditions['methods']))
    
    ### Need to interpret conditions['logic']
    
    # Stub for testing
    return calls[0]
