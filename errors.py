class FunctionError(Exception):
    msg = "Function error."
    status = 422
    
    def __init__(self, info = ""):
        self.info = info
        
class TetherError(Exception):
    msg = "Tether error."
    status = 422
    
    def __init__(self, info = None, request_info = ""):
        self.info = info
        self.request_info = request_info