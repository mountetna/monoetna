from ..list.manipulations import flatten

def parseModelAttr(given: str):
    parts = given.split("#")
    out = {"model": parts[0]}
    if "$" in parts[1]:
        # Table
        att_split = parts[1].split("$")
        out["target_attribute"] = att_split[0]
        out["slice_attribute"] = att_split[1]
        out["slice_content"] = att_split[2]
    else :
        # Non-Table
        out["attribute"]=parts[1]
    return out

def _buildTargetPath_nonTable(targs: str, project_data):
    # Given a model#attribute pair and a project_data definition,
    # Creates the portion of a magma query that can be appended after
    # the [start_model, "::all"]-like initial predicate
    paths = project_data['seq_to_model_paths']
    return flatten( [paths[targs['model']], targs['attribute']] )

def _buildTargetPath_Table(targs: str, project_data):
    # Given a model#attribute pair and a project_data definition,
    # Creates the portion of a magma query that can be appended after
    # the [start_model, "::all"]-like initial predicate
    paths = project_data['seq_to_model_paths']
    return paths[targs['model']] + [
            [targs['slice_attribute'],"::equals",targs['slice_content']],
            "::first", targs["target_attribute"]
        ]

def buildTargetPath(target: str, project_data):
    targs = parseModelAttr(target)
    if len(targs) > 2:
        return _buildTargetPath_Table(targs, project_data)
    return _buildTargetPath_nonTable(targs, project_data)
