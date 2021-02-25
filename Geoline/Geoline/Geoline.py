from typing import List, Dict
from ..Geoline.seqTemplate import *
from functools import partial
from magby import Magby
'''
TODO 
parser of sample name
extractor of md5
directory browser:
    extractor of the file names 
'''

class GeolineError(Exception):
    '''Errors corresponding to misuse of Geoline'''

class Geoline:
    def __init__(self, url: str, token: str) -> None:
        self._url = url
        self._token = token

    def seqWorkflows(self, assay: str) -> None:
        templates = [samplesSection, processedFilesSection, rawFilesSection]
        self.metaDataSections = [self._selectWorkflow(x, assay) for x in templates]


    # Private functions
    def _selectWorkflow(self, template: Callable, assay: str) -> Dict:
        allowedAssays = ['rna_seq', 'dna_seq', 'sc_seq']
        if assay not in allowedAssays:
            raise GeolineError(f'Geoline.selectWorkflow(): unrecognized assay {assay}'
                               f'allowed assays are {allowedAssays}')
        workflow = template(assay)
        updated = self._updater(workflow)
        return updated


    def _updater(self, section: Dict, prev: int=0, curr: int=0) -> Dict:
        if (prev < 0) | (curr >= len(section)):
            return section
        attributeMap = list(section.items())[curr]
        if attributeMap[1] == '':
            return self._updater(section, curr, curr + 1)

        if attributeMap[0] == 'characteristics':
            answer = attributeMap[1](d={})
            return self._updater(dict(section, **{attributeMap[0]: answer}), curr, curr + 1)
        answer = askAttribute(attributeMap[0], attributeMap[1])
        switch = {
            'y': partial(self._updater, prev=curr, curr=curr + 1),
            '0': partial(self._updater, prev=prev - 1, curr=curr - 1),
            'STOP': lambda x: x
        }
        default = lambda x: self._updater(dict(x, **{attributeMap[0]: answer}), prev=curr, curr=curr + 1)
        executor = switch.get(answer, default)
        return executor(section)

    #def _handleValue(self,):










