from archimedes.functions.dataflow import input_json, input_var, output_path
from archimedes.functions.utils import pandas as pd

cluster_meta=input_var('cluster_meta')
discrete_metadata_summary=input_json('discrete_metadata_summary')

pd.DataFrame(
    {
        cluster_meta: discrete_metadata_summary[cluster_meta]
    }
).to_json(output_path('blank_annots.json'))