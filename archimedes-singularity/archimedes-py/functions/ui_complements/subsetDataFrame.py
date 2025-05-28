import pandas as pd
import numpy as np
from ..plotting.utils import _is_logical, _is_continuous, _which_rows
from ..utils import _can_numeric, _can_logical, _all, _to_logical, _can_str
from typing import List
from collections.abc import Iterable

def _check_constraint_format(constraint, index, data_frame):
    issues = []
    type = "unknown"

    if not "col" in constraint:
        issues.append("'col' missing")
    if not "def" in constraint:
        issues.append("'def' missing")
    if  not "logic" in constraint:
        issues.append("'logic' missing")
    elif (index>1 and not constraint['logic'] in ["AND", "OR"]):
        issues.append("'logic' not set to AND or OR")

    # End early
    if (len(issues)>0):
        raise Exception(f'For index-{index+1} constraint: {";".join(issues)}')

    if _is_continuous(data_frame[constraint['col']]):
        # Numeric data target
        type = "numeric"
        if not constraint['def'][0] in ["exactly", "above"]:
            issues.append("numeric 'def' must start with 'exactly' or 'above'")
        if not constraint['def'][2] in ["exactly", "below"] :
            issues.append("numeric 'def' 3rd element must be 'exactly' or 'below'")
        if not _can_numeric(constraint['def'][1]):
            issues.append("numeric 'def' 2nd element must be numeric")
        if not _can_numeric(constraint['def'][3]):
            issues.append("numeric 'def' 4th element must be numeric")
    elif _is_logical(data_frame[constraint['col']]) or _all(_can_logical(data_frame[constraint['col']])):
        # boolean target, actually boolean
        type = "boolean"
        if len(constraint['def'])!= 1:
            issues.append("boolean 'def' must be only one element")
        if not _can_logical(constraint['def'][0], False):
            issues.append("boolean 'def' must equate to [True] or [False]")
    else:
        # string target
        type = "string"
        if len(constraint['def']) < 1 or not _all(_can_str(constraint['col'])):
            issues.append("string 'def' must contain only strings")

    if (len(issues)>0):
        raise Exception(f'For index-{index+1} constraint: {";".join(issues)}')

    return type

def _interpret_constraint(col, _def, type, data_frame):
    target_data = data_frame[col]

    if (type == "numeric"):
        if _def[0]=="above":
            left = target_data > _def[1]
        else:
            left = target_data >= _def[1]
        if _def[2]=="below":
            right = target_data < _def[3]
        else:
            right = target_data <= _def[3]
        out = left and right
    if (type == "boolean"):
        match = _to_logical(_def[0])

        if _is_logical(target_data):
            out = target_data.apply(lambda x: x==match)
        else:
            out = target_data.apply(lambda x: _to_logical(x)==match)
    elif (type == "string"):
        out = target_data.apply(lambda x: x in _def)
    else:
        raise Exception('Unimplemented data-type')

    return out

def subsetDF_index_targets(data_frame, constraints):
    '''
    This function converts from the json-dict returned by the subsetDataFrame UI component,
    into the matching indexes of the target dataframe.
    '''

    raw_calls = []
    logic = []
    for i,v in enumerate(constraints):
        type = _check_constraint_format(v, i, data_frame)
        raw_calls.append(_interpret_constraint(v['col'], v['def'], type, data_frame))
        logic.append(v['logic'])

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
    logic = logic[1:]
    ands = [[0]]
    for i,v in enumerate(logic):
        if v=="AND":
            ands[len(ands)-1].append(i+1)
        elif v=="OR":
            ands.append([i+1])
    and_calls = list(map(lambda this: map_by_2s(combine_ands, [raw_calls[i] for i in this]), ands))
    return _which_rows(map_by_2s(combine_ors, and_calls), data_frame)
