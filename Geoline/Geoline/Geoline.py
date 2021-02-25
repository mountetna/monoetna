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

    def selectWorkflow(self, template: Callable, assay: str):
        allowedAssays = ['rna_seq', 'dna_seq', 'sc_seq']
        if assay not in allowedAssays:
            raise GeolineError(f'Geoline.selectWorkflow(): unrecognized assay {assay}'
                               f'allowed assays are {allowedAssays}')
        allowedTemplates = [samplesSection, processedFilesSection, rawFilesSection]
        if template not in allowedTemplates:
            raise GeolineError(f'Geoline.selectWorkflow(): unrecognized type of workflow {template}'
                               f'supported workflows are {allowedTemplates}')
        workflow = template(assay)
        updated = self._updater(workflow)
        return updated


    def _updater(self, section: dict, prev: int=0, curr: int=0):
        if (prev < 0) | (curr > len(section)):
            return section
        attributeMap = list(section.items())[curr]
        if attributeMap[1] == '':
            return self._updater(section, curr, curr + 1)

        aw = askAttribute(attributeMap[0], attributeMap[1])
        switch = {
            'y': partial(self._updater, curr, curr + 1),
            '0': partial(self._updater, prev - 1, curr - 1),
            'STOP': lambda x: x
        }
        default =  lambda x: self._updater(x.update({attributeMap[0]: aw}), curr, curr + 1)
        executor = switch.get(aw, default)
        return executor(section)










