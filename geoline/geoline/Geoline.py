from typing import List, Dict, Callable, Any
import functools
from functools import partial
from pandas import DataFrame
import re

from magby.Magby import Magby

from .seq_template import samples_section
from .TemplateTree import TemplateTree
from .geo_utils import flatten_nested_dict, ask_attribute

"""
TODO 
parser of sample name
extractor of md5
directory browser:
extractor of the file names 
"""


# Decorator for special treatment of specific columns
def special_columns(function_of_column):
    @functools.wraps(function_of_column)
    def special_columns_wrapper(column: List, *args, **kwargs):
        return [function_of_column(x) for x in column]

    return special_columns_wrapper


class GeolineError(Exception):
    """Errors corresponding to misuse of geoline"""


class Geoline:
    def __init__(self, url: str, token: str, project_name: str) -> None:
        """
        Class to auto populate GEO metadata template
        :param url: str. Magma url
        :param token: str. Private token
        :projectName: str. Name of the project
        """
        self._url = url
        self._token = token
        self._magby_instance = Magby(url, token)
        self._project_name = project_name

    def seq_workflow(self, assay: str, primary_model: str, **kwargs) -> DataFrame:
        """
        Function to trigger interactive auto population of the template.
        :param assay: str. Assay type to focus on. Accepted 'rna_seq', 'dna_seq', 'sc_seq'
        :param primary_model: str. Name of the model that contains most of the data.
        :param kwargs: Additional kwargs to pass to Magma class. Used for passing requests.Session with custom proxies
                        and other attributes
        :return: DataFrame with columns corresponding to the columns in the GEO meta data spreadsheet
        """
        template = samples_section
        attribute_maps = self._select_workflow(template, assay)
        model_groups = self._grouper(attribute_maps)
        query = self._construct_multi_model_query(model_groups, primary_model, **kwargs)
        answer = self._magby_instance.query(self._project_name, query, **kwargs)
        if 'errors' in answer.keys():
            raise GeolineError(f'geoline._workflow(): Malformed magma query {answer["errors"]}')
        flat_answer = self._query_wrapper(answer['answer'], answer['format'])
        answer_df = self._organize_answer_df_columns(DataFrame(flat_answer), attribute_maps)
        return answer_df

    @staticmethod
    def _organize_answer_df_columns(answer_df: DataFrame, attr_map: Dict) -> DataFrame:
        organized_df = DataFrame()
        for geo_field in attr_map:
            if attr_map[geo_field] in answer_df.columns:
                organized_df[geo_field] = answer_df[attr_map[geo_field]]
            else:
                organized_df[geo_field] = ['Not in Magma'] * answer_df.shape[0]
        return organized_df

    def _construct_multi_model_query(self, model_group: Dict, primary_model_name: str, **kwargs) -> List:
        """
        Construct a query for all models at once
        :param model_group: Dict. Structure: {Model : {GEO_meta_data_field : [model_name, attr_name]}}
        :param primary_model_name: str. Name of a model considered to be a starting point for the query
        :param kwargs: Dict. For magby (like session)
        :return: List. A Magby query
        """
        if primary_model_name not in model_group.keys():
            raise GeolineError(f'TemplateTree.multiModelQuery(): {primary_model_name} is not in the list of models '
                               f'{list(model_group.keys())}')
        primary_model = model_group.pop(primary_model_name)
        primary_model_attributes = list(self._extract_attributes(primary_model))
        global_template = self._get_global_template(**kwargs)
        template_tree = TemplateTree(global_template)
        for model in model_group:
            if model != '':
                path_to_model = template_tree.traverse_to_model(primary_model_name, model)
                model_ttr_map = self._extract_attributes(model_group[model])
                #   special treatment of FILE and FILE_COLLECTION attr types
                model_ttr_map = [[x, '::all'] if self._is_file(x) else x for x in model_ttr_map]
                path_to_model.append(model_ttr_map)
                primary_model_attributes.append(path_to_model)
        primary_query = [primary_model_name, '::all', primary_model_attributes]
        return primary_query

    @staticmethod
    def _get_attr_type(global_template: Dict, model: str, attribute: str) -> str:
        return global_template['models'][model]['template']['attributes'][attribute]['attribute_type']

    @staticmethod
    def _is_file(attr_type: str) -> bool:
        return bool(re.search('file', attr_type))

    def _select_workflow(self, template: Callable, assay: str) -> Dict:
        allowed_assays = ['rna_seq', 'dna_seq', 'sc_seq']
        if assay not in allowed_assays:
            raise GeolineError(f'geoline.selectWorkflow(): unrecognized assay {assay}'
                               f'allowed assays are {allowed_assays}')
        workflow = template(assay)
        updated = self._updater(workflow)
        return flatten_nested_dict(updated, sep=' ')

    # TODO refactor
    def _updater(self, section: Dict, prev: int = 0, curr: int = 0) -> Dict:
        # Check out of range
        if (curr < 0) | (curr >= len(section)):
            return section
        attribute_map = list(section.items())[curr]

        # skips empty values in the template
        if attribute_map[1] == '':
            return self._updater(section, curr, curr + 1)

        # special treatment of characteristics group
        if attribute_map[0] == 'characteristics':
            answer = attribute_map[1](d={})
            return self._updater(dict(section, **{attribute_map[0]: answer}), curr, curr + 1)

        answer = ask_attribute(attribute_map[0], attribute_map[1])
        switch = {
            'y': partial(self._updater, prev=curr, curr=curr + 1),
            '0': partial(self._updater, prev=prev - 1, curr=curr - 1),
            'STOP': lambda x: x
        }
        default = lambda x: self._updater(dict(x, **{attribute_map[0]: answer}), prev=curr, curr=curr + 1)
        executor = switch.get(answer, default)
        return executor(section)

    def _grouper(self, section: Dict) -> Dict:
        unique_models = set([self._split_model_map(x)[0] for x in section.values()])
        grouped_models = {x:
                              {y: self._split_model_map(section[y]) for y in section if
                               self._split_model_map(section[y])[0] == x}
                          for x in unique_models}
        return grouped_models

    @staticmethod
    def _split_model_map(attr_map: str) -> List:
        if attr_map == '':
            return ['', '']
        else:
            return attr_map.split(':')

    @staticmethod
    def _extract_attributes(attr_map: Dict) -> List:
        return [x[1] for x in attr_map.values()]

    def _get_global_template(self, **kwargs) -> Dict:
        global_template = self._magby_instance.retrieve(self._project_name, 'all', [], 'all', dataType='json', **kwargs)
        return global_template

    def _walk_answer(self, answer_element: List, answer_format: List) -> Dict:
        # Sometimes in magma some attributes are empty lists. This populates all empty lists in the answer with None
        answer_element = [None,
                          []] if not answer_element else answer_element
        flat_answer = {}
        for elem in range(len(answer_format)):
            if isinstance(answer_format[elem], list):
                flat_answer.update(self._walk_answer(answer_element[elem], answer_format[elem]))
            else:
                if isinstance(answer_element[elem], list):
                    reduced = self._reduce_one_to_many(answer_element)
                    return self._walk_answer(reduced, answer_format)
                flat_answer.update({self._convert_model_attr_to_column(answer_format[elem]): answer_element[elem]})
        return flat_answer

    @staticmethod
    def _convert_model_attr_to_column(model_attr: str) -> str:
        """
        A query to magma returns an answer with "columns" in the format {project::model#attr}.
        Need to convert it to {model:attr} here called geoline Format - to further map to GEO metadata fields
        :param model_attr: str. Magma format
        :return: str. geoline format
        """
        drop_project = model_attr.split("::")[1]
        return drop_project.replace("#", ":")

    def _reduce_one_to_many(self, multi_answer: List) -> List:
        if multi_answer == [None, []]:
            return ["None", "None"]
        reshaped = {x[0]: x[1] for x in multi_answer}
        multi_values = []
        while any([len(x) > 0 for x in reshaped.values()]):
            multi_values.append([self._robust_popper(x) for x in reshaped.values()])
        stringed = [self._any_list_to_string(x) for x in multi_values]
        multi_id = [self._any_list_to_string(list(reshaped.keys())), stringed]
        return multi_id

    def _query_wrapper(self, answer: List, fmt: List) -> List:
        out = [self._walk_answer(answer_element, fmt) for answer_element in answer]
        return out

    @staticmethod
    def _robust_popper(collection: List):
        try:
            return collection.pop(0)
        except IndexError:
            return None

    @staticmethod
    def _any_list_to_string(collection: List[Any], sep: str = ';') -> str:
        return sep.join([str(element) for element in collection])
