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

def characteristics(func: Callable, d: Dict) -> Dict:
    addCharacteristic = addAnother()
    while addCharacteristic in ['y', 'Y']:
        updateInfo = askCharacteristics()
        out.update({updateInfo[0]: updateInfo[1]})

    switch = {
        'Y': characteristics(func, d),
        'y': characteristics(func, d),
        'n': '',
        '0': '',
        'STOP': ''
    }
    updateInfo = switch.get(addCharacteristic)
    out.update({updateInfo[0]: updateInfo[1]})
    return out





from geoUtils import askAttribute
def func(somedict):
    somedict.update({x: askAttribute(x, somedict[x]) for x in somedict})