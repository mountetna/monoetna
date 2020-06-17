import importlib
import methods
import copy
import methods
from errors import FunctionError, TetherError
from traceback import format_exc

def execute(request):
    """
    Executing requests using the methods defined in the methods module.
    
    Request dict format:
    {"func" : <function name that is in the methods module>,
     "requests" : <list of dicts, functions and arguments to run in order. Not used if "func" is not equal to "tether">, 
     "args" : <ordered arguments for func, including data input. For tether, these are the initial input arguments>,
     "kwargs" : <key worded arguments for func>
    }
    
    input: request (dict)
    output: <tether function> OR method.func(*args, **kwargs)
    """
    
    args = request.get("args",[])
    kwargs = request.get("kwargs",{})
        
    func_name = request["func"]
    if func_name == "tether":
        try:
            return tether(request["requests"], args, kwargs)
        except FunctionError as e:
            raise TetherError(request_info = e.info)
        except:
            raise
    else:
        try:
            return methods.FUNCTIONS[func_name](*args, **kwargs)
        except:
            raise FunctionError(info = {"func_name" : func_name, "args" : args, "kwargs" : kwargs})
        
def tether(requests, args, kwargs):
    """
    Tethers multiple requests together, and returns the last return value
    input: requests, args, kwargs
    output: prev_return (last return of the tether list of requests)
    """
    
    prev_return = args
    og_requests = copy.deepcopy(requests)
    
    for request in og_requests:
        if isinstance(prev_return, type(None)):
            pass
        else:
            if type(prev_return) not in set([list, tuple, dict, set]):
                prev_return = [prev_return]
            else:
                prev_return = list(prev_return)
            if "args" not in request:
                request["args"] = prev_return
            else:
                request["args"] = prev_return + request["args"]
            
            if "kwargs" not in request:
                request["kwargs"] = kwargs
            else:
                request["kwargs"].update(prev_return)
         
        prev_return = execute(request)
        
    return prev_return