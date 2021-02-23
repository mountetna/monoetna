from typing import Dict, Callable
from geoUtils import *

def samples() -> Dict:
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
    return flatten(out)


def processedFiles(fileName: str,
                   fileType: str,
                   checksum: str) -> Dict:
    out = {
        'file name': fileName,
        'file type': fileType,
        'file checksum': checksum
    }
    return out


def rawFiles(fileName: str,
             fileType: str,
             checksum: str,
             instrumentModel: str,
             pe: str):
    out = processedFiles(fileName, fileType, checksum)
    out.update({
        'instrument model': instrumentModel,
        'single or paired-end': pe
    })
    return out


def peExperiment(interactiveFunc) -> Dict:
    out = {
        'file name 1': 'ATTR',
        'file name 2': 'ATTR'
    }
    out.update({x: interactiveFunc(x, out[x]) for x in out})
    return out


## INTERACTIVE

def characteristics(addAnother: Callable, d: Dict) -> Dict:
    aw = addAnother
    if aw == 'Y':
        updateInfo = askCharacteristics()
        dc = d.copy()
        dc.update({updateInfo[0]: updateInfo[1]})
        characteristics(addAnother, dc)
    else:
        return d




from geoUtils import askAttribute
def func(somedict):
    somedict.update({x: askAttribute(x, somedict[x]) for x in somedict})