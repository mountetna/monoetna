from typing import List, Dict
import functools
from ..Geoline.seqTemplate import *
from functools import partial
from magby.Magby import Magby
from pandas import DataFrame
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


    def seqWorkflows(self, projectName: str, assay: str, **kwargs) -> None:
        templates = [samplesSection, processedFilesSection, rawFilesSection]
        attributeMaps = [self._selectWorkflow(template, assay) for template in templates]
        modelGroups = [self._grouper(attrMap) for attrMap in attributeMaps]
        magmaExtracts = [self._magmaWrapper(projectName, x, **kwargs) for x in modelGroups]



    # Private functions
    def _selectWorkflow(self, template: Callable, assay: str) -> Dict:
        allowedAssays = ['rna_seq', 'dna_seq', 'sc_seq']
        if assay not in allowedAssays:
            raise GeolineError(f'Geoline.selectWorkflow(): unrecognized assay {assay}'
                               f'allowed assays are {allowedAssays}')
        workflow = template(assay)
        updated = self._updater(workflow)
        return flatten(updated, sep=' ')


    # TODO refactor
    def _updater(self, section: Dict, prev: int=0, curr: int=0) -> Dict:
        # Check out of range
        if (prev < 0) | (curr >= len(section)):
            return section
        attributeMap = list(section.items())[curr]

        # skips empty values in the template
        if attributeMap[1] == '':
            return self._updater(section, curr, curr + 1)

        # special treatment of characteristics group
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


    def _grouper(self, section: Dict) -> Dict:
        uniqueModels = set([self._splitModelMap(x)[0] for x in section.values()])
        groupedModels = {x:
                             {y: self._splitModelMap(section[y]) for y in section if self._splitModelMap(section[y])[0] == x}
                         for x in uniqueModels}
        return groupedModels


    def _splitModelMap(self, attrMap: str) -> List:
        if attrMap == '':
            return ['', '']
        else:
            return attrMap.split(':')

    def _magmaExtractor(self, projectName: str, model: str, attributeMaps: Dict, **kwargs) -> DataFrame:
        magby = Magby(self._url, self._token)
        metadata = magby.retrieve(projectName=projectName,
                                  modelName=model,
                                  recordNames='all',
                                  attributeNames = [x[1] for x in attributeMaps.values()],
                                  dataType='meta',
                                  **kwargs)
        columnMaps = {attributeMaps[x][1]: x for x in attributeMaps}
        return metadata.rename(columns=columnMaps)

    def _magmaWrapper(self, projectName: str, modelGroup: Dict, **kwargs) -> List:
        out = [self._magmaExtractor(projectName, model, modelGroup[model], **kwargs) for model in modelGroup]
        return out



## TODO how to create a magma query for a complex model structure defined by the user??????








