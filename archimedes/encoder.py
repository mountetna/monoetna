import json
from pandas import Series, DataFrame
from numpy import int64

class ArchimedesEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, Series):
            return {
                'vector': [
                    { 'label' : label, 'value' : value }
                    for (label,value) in obj.iteritems()
                ]
            }

        if isinstance(obj, DataFrame):
            return {
                'matrix': {
                    'row_names': obj.index.tolist(),
                    'col_names': [col_name for col_name in obj.columns],
                    'rows': obj.values.tolist()
                }
            }

        if isinstance(obj, int64):
            return obj.item()
        # Let the base class default method raise the TypeError
        return json.JSONEncoder.default(self, obj)

