from numpy import random
import re
import pandas

def do_call(what, kwargs = {}):
    return what(**kwargs)