from ..list.manipulations import flatten

def parseModelAttr(given: str):
    parts = given.split("#")
    return {"model": parts[0], "attribute": parts[1]}

def buildTargetPath(target: str, project_data):
    # Given a model#attribute pair and a project_data definition,
    # Creates the portion of a magma query that can be appended after
    # the [start_model, "::all"]-like initial predicate
    targs = parseModelAttr(project_data[target])
    
    path = (project_data['seq_to_model_paths'])[targs['model']]

    return flatten([path, targs['attribute']])
