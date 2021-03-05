from typing import List, Dict
import functools
from functools import partial
from pandas import DataFrame

from magby.Magby import Magby

from ..Geoline.seqTemplate import *
from ..Geoline.TemplateTree import *

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
    def __init__(self, url: str, token: str, projectName: str) -> None:
        self._url = url
        self._token = token
        self._projectName = projectName


    def seqWorkflows(self, assay: str, **kwargs) -> None:
        templates = [samplesSection, processedFilesSection, rawFilesSection]
        attributeMaps = [self._selectWorkflow(template, assay) for template in templates]
        modelGroups = [self._grouper(attrMap) for attrMap in attributeMaps]



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


    def _extractAttributes(self, attrMap: Dict) -> List:
        return [x[1] for x in attrMap.values()]


    def _getGlobalTemplate(self) -> Dict:
        mb = Magby(self._url, self._token)
        globalTemplate = mb.retrieve(self._projectName, 'all', [], 'all', dataType='json')
        return globalTemplate


    def constructMultiModelQuery(self, modelGroup: Dict, primaryModelName: str) -> List:
        '''
        Construct a query for all models at once
        :param modelGroup: Dict. Structure: {Model : {GEO_meta_data_field : [model_name, attr_name]}}
        :param primaryModelName: str. Name of a model considered to be a starting point for the query
        :return: List. A Magby query
        '''
        if primaryModelName not in modelGroup.keys():
            raise GeolineError(f'TemplateTree.multiModelQuery(): {primaryModelName} is not in the list of models '
                               f'{list(modelGroup.keys())}')
        primaryModel = modelGroup.pop(primaryModelName)
        primaryModelAttributes = list(self._extractAttributes(primaryModel))
        templateTree = TemplateTree(self._getGlobalTemplate())
        for model in modelGroup:
            if model != '':
                pathToModel = templateTree.traverseToModel(primaryModelName, model)+['::all']
                modelAttrMap = self._extractAttributes(modelGroup[model])
                pathToModel.append(modelAttrMap)
                primaryModelAttributes.append(pathToModel)
        primaryQuery = [primaryModelName, '::all']
        primaryQuery.append(primaryModelAttributes)
        return primaryQuery


    def _walkAnswer(self, answerElement: List, answerFormat: List, flattenedAnswer: Dict) -> Dict:
        for elementNum in range(len(answerFormat)):
            if isinstance(answerFormat[elementNum], list):
                return self._walkAnswer(answerElement[elementNum], answerFormat[elementNum], flattenedAnswer)
            if isinstance(answerElement[elementNum], list):
                reshapedAnswerElement = [list(x) for x in zip(*answerElement[elementNum])]
                flattenedAnswer.update({answerFormat[elementNum]: reshapedAnswerElement[0]})
                return self._walkAnswer(reshapedAnswerElement[1:], answerFormat[1:], flattenedAnswer)
        flattenedAnswer.update({answerElement[0]: answerFormat[0]})
        return flattenedAnswer




    def _magmaWrapper(self, projectName: str, modelGroup: Dict, **kwargs) -> List:
        out = [self._magmaExtractor(projectName, model, modelGroup[model], **kwargs) for model in modelGroup]
        return out



