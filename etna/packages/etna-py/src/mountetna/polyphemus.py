import dataclasses
import typing
from io import StringIO
from functools import partial
from inspect import isgenerator
from typing import Dict, Optional, List
from serde import serialize, deserialize
from serde.json import from_json, to_json

import pandas as pd

from .etna_base import EtnaClientBase
from requests import HTTPError
from .utils.iterables import batch_iterable
from .metis import Upload, File


class Polyphemus(EtnaClientBase):
    pass
