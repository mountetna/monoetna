def _default_to_if_make_and_logic(this, default, logic = True):
    if (this == "make" and logic):
        return default
    return this
