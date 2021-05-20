from archimedes.functions.dataflow import curl_data, input_json, input_var, output_json
from DataIO import RemoteData

if __name__ == '__main__':
    url = str(input_var('url'))
    h5io = RemoteData(url)
    h5io.data_from_metis_to_memory()
    h5io.write_h5ad('processed_data.h5ad')

    covariates = h5io.extract_covariate_json()
    output_json(covariates, 'covariates.json')