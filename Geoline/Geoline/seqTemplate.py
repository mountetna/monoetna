from typing import Dict, Callable
from ..Geoline.geoUtils import *


def samplesSection(title: str,
                   sourceName: str,
                   organism: str,
                   characteristics: Dict,
                   molecule: str,
                   description: str,
                   processedDataFile: str,
                   rawFile: Dict) -> Dict:
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
    return flatten(out, sep=' ')


def processedFilesSection(fileName: str,
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
    out = processedFilesSection(fileName, fileType, checksum)
    out.update({
        'instrument model': instrumentModel,
        'single or paired-end': pe
    })
    return out


def peExperimentSection(interactiveFunc) -> Dict:
    out = {
        'file name 1': 'ATTR',
        'file name 2': 'ATTR'
    }
    out.update({x: interactiveFunc(x, out[x]) for x in out})
    return out


def characteristics(addAnother: Callable, d: Dict) -> Dict:
    aw = addAnother()
    if aw == 'y':
        updateInfo = askCharacteristics()
        dc = d.copy()
        dc.update({updateInfo[0]: updateInfo[1]})
        return characteristics(addAnother, dc)
    else:
        return d

