from archimedes.functions.dataflow import curl_data, input_path, output_path
from archimedes.functions.scanpy import scanpy as sc

from typing import Callable
from io import BytesIO
import json

class RemoteData(object):
    def __init__(self, url: str) -> None:
        self._url = url
        self._buffer = BytesIO()
        self.anndata = None

    def data_from_metis_to_memory(self) -> None:
        self._buffer.write(curl_data(self._url).content)
        self._buffer.seek(0)
        self.anndata = sc.read_h5ad(self._buffer)
        self._buffer.close()

    def write_h5ad(self, file_name: str) -> None:
        self.anndata.write(file_name)

    def extract_covariate_json(self) -> str:
        # TODO handle other dtypes (not string or categorical) correctly, eg provide thresholds for numerical
        covariates = {
            x: list(self.anndata.obs[x].unique()) for x in self.anndata.obs
        }
        return json.dumps(covariates)


class LocalData(object):
    def __init__(self, data_location: str) -> None:
        self._dir = data_location

    def io_wrapper(self, inner_func: Callable, output_name: str, **kwargs) -> None:
        """
        Higher order function to handle AnnData io and perform transformations
        :param inner_func: Callable. Function to run on AnnData, must return transformed AnnData
        :param kwargs: dict. For inner_func
        :return: None
        """
        scdata = sc.read(input_path(self._dir))
        scdata = inner_func(scdata, **kwargs)
        scdata.write(output_path(output_name))




