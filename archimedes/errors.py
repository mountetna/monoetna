class ArchimedesError(Exception):
    msg = "Archimedes error."
    status = 422
    
    def __init__(self, info = ""):
        self.info = info
