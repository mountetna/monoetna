from functools import partial

from geoline.geoUtils import *


def samplesSection(assay: str) -> Dict:
    out = {
        'title': f'{assay}:tube_name',
        'source name': f'{assay}:biospecimen',
        'organism': 'subject:name',
        'characteristics': partial(characteristics, addAnother=addAnother),
        'molecule': '',
        'description': f'{assay}:notes',
        'processed data file': f'{assay}:gene_expression',
        'raw': f'{assay}:raw_fastqs'
    }
    return out



