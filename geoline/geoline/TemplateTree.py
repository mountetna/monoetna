from typing import List, Dict


class TemplateTree(object):
    def __init__(self, global_template: Dict):
        self.global_template = global_template

    def _ascend_tree(self, travel_path: List) -> List:
        """
        Return ascending path from starting node to a root node
        :param travel_path: List. List with a single element, which is a starting node
        :return: List of nodes
        """
        try:
            parent = self.global_template['models'][travel_path[-1]]['template']['parent']
            return self._ascend_tree(travel_path + [parent])
        except KeyError:
            return travel_path

    def _common_root(self, primary_tree: List, secondary_tree: List) -> str:
        if primary_tree == []:
            return ''
        if primary_tree[0] in secondary_tree:
            return primary_tree[0]
        return self._common_root(primary_tree[1:], secondary_tree)

    def traverse_to_model(self, primary_model: str, secondary_model: str) -> List:
        primary_tree, secondary_tree = [self._ascend_tree([x]) for x in [primary_model, secondary_model]]
        common_root = self._common_root(primary_tree, secondary_tree)
        new_path = primary_tree[1:primary_tree.index(common_root)] + secondary_tree[:(secondary_tree.index(common_root)) + 1][
                                                                 ::-1]
        # Branching models
        if not all(x in primary_tree for x in secondary_tree):
            new_path.append('::all')
        return new_path
