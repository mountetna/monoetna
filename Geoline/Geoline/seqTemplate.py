from ..Geoline.geoUtils import *


def samplesSection(assay: str) -> Dict:
    out = {
        'title': f'{assay}:tube_name',
        'source name': f'{assay}:biospecimen',
        'organism': 'subject:name',
        'characteristics': characteristics(addAnother, {}),
        'molecule': '',
        'description': f'{assay}:notes',
        'processed data file': f'{assay}:gene_expression',
        'raw': ''
    }
    return flatten(out, sep=' ')


def processedFilesSection(assay: str) -> Dict:
    out = {
        'file name': f'{assay}:gene_expression',
        'file type': '',
        'file checksum': ''
    }
    return out


def rawFilesSection(assay: str):
    out = processedFilesSection(assay)
    out.update({
        'instrument model': f'{assay}:sequencer',
        'single or paired-end': ''
    })
    return out

