from functools import partial
from typing import Dict

from .geoUtils import characteristics, addAnother


def samples_section(assay: str) -> Dict:
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



