def get_modules():
    from os.path import dirname, basename, isfile, join
    import glob
    import importlib
    
    return [importlib.import_module("methods.%s" % basename(f)[:-3])
            for f in glob.glob(join(dirname(__file__), "*.py"))
            if isfile(f) and basename(f) != "__init__.py"]

def get_functions():
    from os.path import dirname, basename, isfile, join
    import glob
    import importlib
  
    return dict([func_key_value 
            for module in get_modules()
            for func_key_value in [[f, getattr(module, f)] 
                                    for f in [m for m in dir(module)
                                              if (getattr(module, m).__name__ is "pythia_func" 
                                                  if "__name__" in dir(getattr(module, m)) else False)]]])
    
#MODULES = get_modules()
FUNCTIONS = get_functions()