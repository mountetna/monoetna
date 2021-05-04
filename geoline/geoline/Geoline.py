from typing import List, Dict, Callable
import functools
from functools import partial
from pandas import DataFrame
import re

from magby.Magby import Magby

from .seqTemplate import samples_section
from .TemplateTree import TemplateTree
from .geoUtils import flatten, askAttribute

'''
TODO 
parser of sample name
extractor of md5
directory browser:
    extractor of the file names 
'''


# Decorator for special treatment of specific columns
def special_columns(function_of_column):
    @functools.wraps(function_of_column)
    def special_columns_wrapper(column: List, *args, **kwargs):
        return [function_of_column(x) for x in column]

    return special_columns_wrapper


class GeolineError(Exception):
    '''Errors corresponding to misuse of geoline'''


class Geoline:
    def __init__(self, url: str, token: str, project_name: str) -> None:
        '''
        Class to auto populate GEO metadata template
        :param url: str. Magma url
        :param token: str. Private token
        :projectName: str. Name of the project
        '''
        self._url = url
        self._token = token
        self._magbyInstance = Magby(url, token)
        self._projectName = project_name

    def seqWorkflow(self, assay: str, primary_model: str, **kwargs) -> DataFrame:
        '''
        Function to trigger interactive auto population of the template.
        :param assay: str. Assay type to focus on. Accepted 'rna_seq', 'dna_seq', 'sc_seq'
        :param primary_model: str. Name of the model that contains most of the data.
        :param kwargs: Additional kwargs to pass to Magma class. Used for passing requests.Session with custom proxies
                        and other attributes
        :return: DataFrame with columns corresponding to the columns in the GEO meta data spreadsheet
        '''
        template = samples_section
        attribute_maps = self._select_workflow(template, assay)
        modelGroups = self._grouper(attribute_maps)
        query = self._construct_multi_model_query(modelGroups, primary_model, **kwargs)
        answer = self._magbyInstance.query(self._projectName, query, **kwargs)
        if 'errors' in answer.keys():
            raise GeolineError(f'geoline._workflow(): Malformed magma query {answer["errors"]}')
        flatAnswer = self._query_wrapper(answer['answer'], answer['format'])
        answer_df = self._organize_answer_df_columns(DataFrame(flatAnswer), attribute_maps)
        return answer_df

    def _organize_answer_df_columns(self, answer_df: DataFrame, attr_map: Dict) -> DataFrame:
        organizedDF = DataFrame()
        for geoField in attr_map:
            if attr_map[geoField] in answer_df.columns:
                organizedDF[geoField] = answer_df[attr_map[geoField]]
            else:
                organizedDF[geoField] = ['Not in Magma'] * answer_df.shape[0]
        return organizedDF

    def _construct_multi_model_query(self, model_group: Dict, primary_model_name: str, **kwargs) -> List:
        '''
        Construct a query for all models at once
        :param model_group: Dict. Structure: {Model : {GEO_meta_data_field : [model_name, attr_name]}}
        :param primary_model_name: str. Name of a model considered to be a starting point for the query
        :param kwargs: Dict. For magby (like session)
        :return: List. A Magby query
        '''
        if primary_model_name not in model_group.keys():
            raise GeolineError(f'TemplateTree.multiModelQuery(): {primary_model_name} is not in the list of models '
                               f'{list(model_group.keys())}')
        primaryModel = model_group.pop(primary_model_name)
        primaryModelAttributes = list(self._extract_attributes(primaryModel))
        globalTemplate = self._get_global_template(**kwargs)
        templateTree = TemplateTree(globalTemplate)
        for model in model_group:
            if model != '':
                pathToModel = templateTree.traverseToModel(primary_model_name, model)
                modelAttrMap = self._extract_attributes(model_group[model])
                #   special treatment of FILE and FILE_COLLECTION attr types
                modelAttrMap = [[x, '::all'] if self._is_file(x) else x for x in modelAttrMap]
                pathToModel.append(modelAttrMap)
                primaryModelAttributes.append(pathToModel)
        primaryQuery = [primary_model_name, '::all']
        primaryQuery.append(primaryModelAttributes)
        return primaryQuery

    def _get_attr_type(self, global_template: Dict, model: str, attribute: str) -> str:
        return global_template['models'][model]['template']['attributes'][attribute]['attribute_type']

    def _is_file(self, attr_type: str) -> bool:
        return bool(re.search('file', attr_type))

    def _select_workflow(self, template: Callable, assay: str) -> Dict:
        allowedAssays = ['rna_seq', 'dna_seq', 'sc_seq']
        if assay not in allowedAssays:
            raise GeolineError(f'geoline.selectWorkflow(): unrecognized assay {assay}'
                               f'allowed assays are {allowedAssays}')
        workflow = template(assay)
        updated = self._updater(workflow)
        return flatten(updated, sep=' ')

    # TODO refactor
    def _updater(self, section: Dict, prev: int = 0, curr: int = 0) -> Dict:
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
        uniqueModels = set([self._split_model_map(x)[0] for x in section.values()])
        groupedModels = {x:
                             {y: self._split_model_map(section[y]) for y in section if
                              self._split_model_map(section[y])[0] == x}
                         for x in uniqueModels}
        return groupedModels

    def _split_model_map(self, attr_map: str) -> List:
        if attr_map == '':
            return ['', '']
        else:
            return attr_map.split(':')

    def _extract_attributes(self, attr_map: Dict) -> List:
        return [x[1] for x in attr_map.values()]

    def _get_global_template(self, **kwargs) -> Dict:
        globalTemplate = self._magbyInstance.retrieve(self._projectName, 'all', [], 'all', dataType='json', **kwargs)
        return globalTemplate

    def _walkAnswer(self, answer_element: List, answer_format: List) -> Dict:
        # Sometimes in magma some attributes are empty lists. This populates all empty lists in the answer with None
        answer_element = [None,
                          []] if not answer_element else answer_element
        flatAnswer = {}
        for elem in range(len(answer_format)):
            if isinstance(answer_format[elem], list):
                flatAnswer.update(self._walkAnswer(answer_element[elem], answer_format[elem]))
            else:
                if isinstance(answer_element[elem], list):
                    reduced = self._reduce_one_to_many(answer_element)
                    return self._walkAnswer(reduced, answer_format)
                flatAnswer.update({self._convert_model_attr_to_column(answer_format[elem]): answer_element[elem]})
        return flatAnswer

    def _convert_model_attr_to_column(self, model_attr: str) -> str:
        '''
        A query to magma returns an answer with "columns" in the format {project::model#attr}.
        Need to convert it to {model:attr} here called geoline Format - to further map to GEO metadata fields
        :param model_attr: str. Magma format
        :return: str. geoline format
        '''
        dropProject = model_attr.split("::")[1]
        return dropProject.replace("#", ":")

    def _reduce_one_to_many(self, multi_answer: List) -> List:
        reshaped = [list(x) for x in zip(*multi_answer)]
        stringed = ['; '.join([str(word) for word in x]) for x in reshaped]
        return stringed

    def _query_wrapper(self, answer: List, format: List) -> List:
        out = [self._walkAnswer(answerElement, format) for answerElement in answer]
        return out
