from archimedes.functions.dataflow import input_json, input_var, input_bool, input_path, output_path
from archimedes.functions.scanpy import scanpy as sc
from archimedes.functions.utils import pandas as pd
from archimedes.functions.list import unique

## Read in data
DEdf = pd.read_csv(input_path('full_diffexp.csv'))

min_abs_fc = input_var('min_abs_fc')
max_pval = input_var('max_pval')
min_pct = input_var('min_pct')
pos_only = input_bool('pos_only')

# DEdf_filtered = 

## Output Filtered Results
DEdf.to_csv(output_path('filtered_diffexp.csv'))