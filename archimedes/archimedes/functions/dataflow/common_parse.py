def parseModelAttr(given: str):
    parts = given.split(":")
    return {"model": parts[0], "attribute": parts[1]}