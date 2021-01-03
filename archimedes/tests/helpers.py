from request_handler import resolve
from encoder import ArchimedesEncoder
import json

def resolve_json(manifest):
    payload = resolve(manifest)

    return json.loads(json.dumps(payload, cls=ArchimedesEncoder))

def labels(v):
    return [ k['label'] for k in v ]

def values(v):
    return [ k['value'] for k in v ]
