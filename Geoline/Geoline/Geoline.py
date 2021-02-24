from typing import List, Dict
from ..Geoline.workflows import *

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

    def selectWorkflow(self, assay: str):
        allowedAssays = ['rna_seq', 'dna_seq', 'sc_seq']
        samples = sampleMapAttr(assay)
        for field, magmaAttr in samples.items():
            aw = askAttribute(field, magmaAttr)
            if aw == 'Y':
                continue
            else:
                samples.update({field: aw})

    # THIS
    def _updater(self, prev, curr):
        pass









