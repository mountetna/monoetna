from typing import List, Dict
import functools
from functools import partial
from pandas import DataFrame
import re

from magby.Magby import Magby

from geoline.seqTemplate import *
from geoline.TemplateTree import *

'''
TODO 
parser of sample name
extractor of md5
directory browser:
    extractor of the file names 
'''

# Decorator for special treatment of specific columns
def specialColumns(functionOfColumn):
    @functools.wraps(functionOfColumn)
    def specialColumnsWrapper(column: List, *args, **kwargs):
        return [functionOfColumn(x) for x in column]
    return specialColumnsWrapper

class GeolineError(Exception):
    '''Errors corresponding to misuse of geoline'''


class Geoline:
    def __init__(self, url: str, token: str, projectName: str) -> None:
        '''
        Class to auto populate GEO metadata template
        :param url: str. Magma url
        :param token: str. Private token
        :projectName: str. Name of the project
        '''
        self._url = url
        self._token = token
        self._magbyInstance = Magby(url, token)
        self._projectName = projectName


    def seqWorkflow(self, assay: str, primaryModel: str, **kwargs) -> DataFrame:
        '''
        Function to trigger interactive auto population of the template.
        :param assay: str. Assay type to focus on. Accepted 'rna_seq', 'dna_seq', 'sc_seq'
        :param primaryModel: str. Name of the model that contains most of the data.
        :param kwargs: Additional kwargs to pass to Magma class. Used for passing requests.Session with custom proxies
                        and other attributes
        :return: DataFrame with columns corresponding to the columns in the GEO meta data spreadsheet
        '''
        template = samplesSection
        attributeMaps = self._selectWorkflow(template, assay)
        modelGroups = self._grouper(attributeMaps)
        query = self._constructMultiModelQuery(modelGroups, primaryModel, **kwargs)
        answer = self._magbyInstance.query(self._projectName, query, **kwargs)
        if 'errors' in answer.keys():
            raise GeolineError(f'geoline._workflow(): Malformed magma query {answer["errors"]}')
        flatAnswer = self._queryWrapper(answer['answer'], answer['format'])
        answerDF = self._organizeAnswerDFColumns(DataFrame(flatAnswer), attributeMaps)
        return answerDF


    def _organizeAnswerDFColumns(self, answerDF: DataFrame, attrMap: Dict) -> DataFrame:
        organizedDF = DataFrame()
        for geoField in attrMap:
            if attrMap[geoField] in answerDF.columns:
                organizedDF[geoField] = answerDF[attrMap[geoField]]
            else:
                organizedDF[geoField] = ['Not in Magma'] * answerDF.shape[0]
        return organizedDF


    def _constructMultiModelQuery(self, modelGroup: Dict, primaryModelName: str, **kwargs) -> List:
        '''
        Construct a query for all models at once
        :param modelGroup: Dict. Structure: {Model : {GEO_meta_data_field : [model_name, attr_name]}}
        :param primaryModelName: str. Name of a model considered to be a starting point for the query
        :param kwargs: Dict. For magby (like session)
        :return: List. A Magby query
        '''
        if primaryModelName not in modelGroup.keys():
            raise GeolineError(f'TemplateTree.multiModelQuery(): {primaryModelName} is not in the list of models '
                               f'{list(modelGroup.keys())}')
        primaryModel = modelGroup.pop(primaryModelName)
        primaryModelAttributes = list(self._extractAttributes(primaryModel))
        globalTemplate = self._getGlobalTemplate(**kwargs)
        templateTree = TemplateTree(globalTemplate)
        for model in modelGroup:
            if model != '':
                pathToModel = templateTree.traverseToModel(primaryModelName, model)
                modelAttrMap = self._extractAttributes(modelGroup[model])
                #   special treatment of FILE and FILE_COLLECTION attr types
                modelAttrMap = [[x, '::all'] if self._isFile(x) else x for x in modelAttrMap]
                pathToModel.append(modelAttrMap)
                primaryModelAttributes.append(pathToModel)
        primaryQuery = [primaryModelName, '::all']
        primaryQuery.append(primaryModelAttributes)
        return primaryQuery


    def _getAttrType(self, globalTemplate: Dict, model: str, attribute: str) -> str:
        return globalTemplate['models'][model]['template']['attributes'][attribute]['attribute_type']

    def _isFile(self, attrType: str) -> bool:
        return bool(re.search('file', attrType))

    def _selectWorkflow(self, template: Callable, assay: str) -> Dict:
        allowedAssays = ['rna_seq', 'dna_seq', 'sc_seq']
        if assay not in allowedAssays:
            raise GeolineError(f'geoline.selectWorkflow(): unrecognized assay {assay}'
                               f'allowed assays are {allowedAssays}')
        workflow = template(assay)
        updated = self._updater(workflow)
        return flatten(updated, sep=' ')


    # TODO refactor
    def _updater(self, section: Dict, prev: int=0, curr: int=0) -> Dict:
        # Check out of range
        if (curr < 0) | (curr >= len(section)):
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


    def _getGlobalTemplate(self, **kwargs) -> Dict:
        globalTemplate = self._magbyInstance.retrieve(self._projectName, 'all', [], 'all', dataType='json', **kwargs)
        return globalTemplate


    def _walkAnswer(self, answerElement: List, answerFormat: List) -> Dict:
        answerElement = [None, []] if not answerElement else answerElement # Sometimes in magma some attributes are empty lists. This populates all empty lists in the answer with None
        flatAnswer = {}
        for elem in range(len(answerFormat)):
            if isinstance(answerFormat[elem], list):
                flatAnswer.update(self._walkAnswer(answerElement[elem], answerFormat[elem]))
            else:
                if isinstance(answerElement[elem], list):
                    reduced = self._reduceOneToMany(answerElement)
                    return self._walkAnswer(reduced, answerFormat)
                flatAnswer.update({self._convertModelAttrToColumn(answerFormat[elem]): answerElement[elem]})
        return flatAnswer


    def _convertModelAttrToColumn(self, modelAttr: str) -> str:
        '''
        A query to magma returns an answer with "columns" in the format {project::model#attr}.
        Need to convert it to {model:attr} here called geoline Format - to further map to GEO metadata fields
        :param modelAttr: str. Magma format
        :return: str. geoline format
        '''
        dropProject = modelAttr.split("::")[1]
        return dropProject.replace("#", ":")



    def _reduceOneToMany(self, multiAnswer: List) -> List:
        reshaped = [list(x) for x in zip(*multiAnswer)]
        stringed = ['; '.join([str(word) for word in x]) for x in reshaped]
        return stringed


    def _queryWrapper(self, answer: List, format: List) -> List:
        out = [self._walkAnswer(answerElement, format) for answerElement in answer]
        return out



