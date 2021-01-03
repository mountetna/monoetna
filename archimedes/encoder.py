import json
from pandas import Series;

class ArchimedesEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, Series):
            return [
                { 'label' : label, 'value' : value }
                for (label,value) in obj.iteritems()
            ]
        # Let the base class default method raise the TypeError
        return json.JSONEncoder.default(self, obj)

