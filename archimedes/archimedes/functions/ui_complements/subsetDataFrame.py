import pandas as pd
import numpy as np
from archimedes.functions.plotting.utils import _is_logical, _is_continuous, _which_rows
from typing import List

def subsetDF_index_targets(data_frame, conditions):
    '''
    This function converts from the json-dict returned by the subsetDataFrame UI component,
    into the matching indexes of the target dataframe.
    '''
    
    if "methods" not in conditions.keys() or "logic" not in conditions.keys():
        raise Exception("Subsetting conditions are not formatted properly.")
    
    # Parse when each condition is True/False.
    def interpret_method(method: List[str]):
        if method[0] is None or method[0]==np.NaN:
            raise Exception("A subsetting condition is incomplete.")
        target_data = data_frame[method[0]]
        # Numeric data
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
        # Logical data
        if _is_logical(target_data):
            if len(method)!=2 or not method[1].lower() in ["true", "false"]:
                raise Exception("A subsetting condition for logical-type data was left incomplete.")
            match = False
            if method[1].lower()=="true":
                match = True
            return [a==match for a in target_data]
        # Remainder = String data
        if len(method)<2:
            raise Exception("A subsetting condition for string-type data was left incomplete, or a condition targetting an unimplemented data type was left in place.")
        return [a in method[1:] for a in target_data]
    raw_calls = list(map(interpret_method, conditions['methods']))
    
    # Combine all calls by and/or logic
    def combine_ors(clause1, clause2):
        return [any([a,b]) for a,b in zip(clause1, clause2)]
    def combine_ands(clause1, clause2):
        return [all([a,b]) for a,b in zip(clause1, clause2)]
    def map_by_2s(fxn, calls: List[List[bool]]):
        out = calls[0]
        i=1
        while i < len(calls):
            out = fxn(out, calls[i])
            i += 1
        return out        
    logic = conditions['logic']
    ands = [[0]]
    if logic[0]!=[]:
        for i in range(len(logic)):
            if logic[i]==["and"]:
                ands[len(ands)-1].append(i+1)
            elif logic[i]==["or"]:
                ands.append([i+1])
    and_calls = list(map(lambda this: map_by_2s(combine_ands, [raw_calls[i] for i in this]), ands))
    return _which_rows(map_by_2s(combine_ors, and_calls), data_frame)
