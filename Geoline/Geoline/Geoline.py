from typing import List, Dict
from ..Geoline.workflows import *
from functools import partial
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
        updated = self._updater(samples)


    ### TEST THIS STUFF
    def _updater(self, section: dict, prev: int=0, curr: int=0):
        if (prev < 0) | (curr > len(section)):
            return section
        curr = section.items()[curr]
        aw = askAttribute(curr[0], curr[1])
        switch = {
            'Y': partial(self._updater, curr, curr + 1),
            '0': partial(self._updater, prev - 1, curr - 1),
            'STOP': lambda x: x
        }
        default =  lambda x: self._updater(x.update({curr[0]: curr[1]}), curr, curr + 1)
        executor = switch.get(aw, default)
        return executor(section)










