from archimedes.functions.dataflow import input_json, input_var, input_bool, input_path, output_path
from archimedes.functions.scanpy import scanpy as sc
from archimedes.functions.utils import pandas as pd
from archimedes.functions.list import unique

## Read in data
DEdf = pd.read_csv(input_path('full_diffexp.csv'), header=0, index_col=0)

# LogFC & pos only
DEdf_filt = DEdf
pos_only = input_bool('pos_only')
if pos_only:
    DEdf_filt = DEdf[ DEdf['logfoldchanges'] > 0 ]
min_abs_fc = float(input_var('min_abs_fc'))
DEdf_filt = DEdf_filt[ abs(DEdf_filt['logfoldchanges']) >= min_abs_fc ]

# pval_adj
max_pval = float(input_var('max_pval'))
DEdf_filt = DEdf_filt[ DEdf_filt['pvals_adj'] <= max_pval ]

# min_pct
min_pct = float(input_var('min_pct'))
pct_cols = [s for s in DEdf.keys() if s.startswith("pts_")]
# Algorithm choice: OR filter, a.k.a. either group meets the cutoff
idx_passed = [
    i for i in list(DEdf_filt.index) if (
    list(DEdf[pct_cols[0]])[i] >= min_pct or list(DEdf[pct_cols[1]])[i] >= min_pct
    )]
DEdf_filt = DEdf_filt.loc[idx_passed,:]

## Output Filtered Results
DEdf_filt.to_csv(output_path('filtered_diffexp.csv'))