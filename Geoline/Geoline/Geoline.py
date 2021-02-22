from typing import List, Dict


from magby import Magby
'''
call in the CLI (url, token)
prompt select project
autopopulate "Series"
prompt assays or cancel
for each column in "Sample" suggest model attr
    confirm, suggest your own, back, cancel
"Processed data files" 
    confirm, suggest your own, back, cancel
"RAW data files" 
    confirm, suggest your own, back, cancel
"Pair End" 
    confirm, suggest your own, back, cancel
    
START
'''

class Geoline:
    def __init__(self, url: str, token: str) -> None:
        self._url = url
        self._token = token


    def _samplesSection(self,
                        title: str,
                        sourceName: str,
                        organism: str,
                        characteristics: Dict,
                        molecule: str,
                        description: str,
                        processedDataFile: str,
                        rawFile: Dict):
        out = {
            'title': title,
            'source name': sourceName,
            'organism': organism,
            'characteristics': characteristics,
            'molecule': molecule,
            'description': description,
            'processed data file': processedDataFile,
            'raw file': rawFile
        }
        return self._flatten(out)


    def _flatten(self, dictionary, parent_key='', sep=' '):
        items = []
        for currKey, currVal in dictionary.items():
            new_key = parent_key + sep + currKey if parent_key else currKey
            if isinstance(currVal, dict):
                items.extend(self._flatten(currVal, new_key, sep=sep).items())
            else:
                items.append((new_key, currVal))
        return dict(items)





