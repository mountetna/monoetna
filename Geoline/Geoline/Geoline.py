from typing import List, Dict


from magby import Magby
'''
call in the CLI (url, token)
prompt select project
autopopulate "Series"
prompt assays or cancel
for each column in "Sample" suggest model attr
    confirm, suggest your own, back, cancel
"Processed data files" 
    confirm, suggest your own, back, cancel
"RAW data files" 
    confirm, suggest your own, back, cancel
"Pair End" 
    confirm, suggest your own, back, cancel
    
START
'''

class Geoline:
    def __init__(self, url: str, token: str) -> None:
        self._url = url
        self._token = token



