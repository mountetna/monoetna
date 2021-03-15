from typing import List, Dict
import functools
from functools import partial

from magby.Magby import Magby

from geoline.geoline.seqTemplate import *
from geoline.geoline.TemplateTree import *
import geoline.geoline.asperaSDK.transfer_pb2 as transfer_manager
import geoline.geoline.asperaSDK.transfer_pb2_grpc as transfer_manager_grpc


'''
From requests:

then this stream to the SDK

resume upload: aspera if file exist return it's size in bytes, then read from offset and 
'''


class Geospera(object):
    def __init__(self, url: str, token: str) -> None:
        self._url = url
        self._token = token


    def upload(self, fileName):
        pass
