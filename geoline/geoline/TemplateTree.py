from typing import List, Dict

class TemplateTree(object):
    def __init__(self, globalTemplate: Dict):
        self.globalTemplate = globalTemplate


    def _ascendTree(self, travelPath: List) -> List:
        '''
        Return ascending path from starting node to a root node
        :param travelPath: List. List with a single element, which is a starting node
        :return: List of nodes
        '''
        try:
            parent = self.globalTemplate['models'][travelPath[-1]]['template']['parent']
            return self._ascendTree(travelPath + [parent])
        except KeyError:
            return travelPath


    def _commonRoot(self, primaryTree: List, secondaryTree: List) -> str:
        if primaryTree == []:
            return ''
        if primaryTree[0] in secondaryTree:
            return primaryTree[0]
        return self._commonRoot(primaryTree[1:], secondaryTree)


    def traverseToModel(self, primaryModel: str, secondaryModel: str) -> List:
        primaryTree, secondaryTree = [self._ascendTree([x]) for x in [primaryModel, secondaryModel]]
        commonRoot = self._commonRoot(primaryTree, secondaryTree)
        newPath = primaryTree[1:primaryTree.index(commonRoot)] + secondaryTree[:(secondaryTree.index(commonRoot))+1][::-1]
        # Branching models
        if not all(x in primaryTree for x in secondaryTree):
            newPath.append('::all')
        return newPath












