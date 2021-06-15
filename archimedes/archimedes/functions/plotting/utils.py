def default_to_if_None_and_logic(this, fill, logic = True):
    if (this == None and logic):
        return fill
    return this
