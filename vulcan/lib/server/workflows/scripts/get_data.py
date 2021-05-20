from archimedes.functions.dataflow import curl_data, input_json, input_var
from archimedes.functions.scanpy import scanpy as sc
from io import BytesIO


class ProcessedData(object):
    def __init__(self, url: str) -> None:
        self._url = url
        self._buffer = BytesIO()
        self.anndata = None

    def get_data(self) -> None:
        self._buffer.write(curl_data(self._url).content)
        self._buffer.seek(0)
        self.anndata = sc.read_h5ad(self._buffer)
        self._buffer.close()

    def write_data(self, file_name: str):
        self.anndata.write(file_name)


if __name__ == '__main__':
    url = str(input_var('url'))
    processed_data = ProcessedData(url)
    processed_data.get_data()
    processed_data.write_data('processed_data.h5ad')