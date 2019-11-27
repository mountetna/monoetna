from methods import algebra
from methods import analysis
import importlib

def get_func(func_name):
    if len(func_name.split(".")) == 2:
        [module, func] = func_name.split(".")
    else:
        [module, func] = [func_name, "func"]
    
    # methods.module should be already loaded above...
    try:
        return getattr(importlib.import_module("%s.%s"%("methods", module)), func)