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
    DEdf_filt = DEdf_filt[ DEdf_filt['logfoldchanges'] > 0 ]
min_abs_fc = float(input_var('min_abs_fc'))
DEdf_filt = DEdf_filt[ abs(DEdf_filt['logfoldchanges']) >= min_abs_fc ]

# pval_adj
max_pval = float(input_var('max_pval'))
DEdf_filt = DEdf_filt[ DEdf_filt['pvals_adj'] <= max_pval ]

# min_pct
min_pct = float(input_var('min_pct'))
# Algorithm choice: OR filter, a.k.a. either group meets the cutoff
pct_cols = [s for s in DEdf_filt.keys() if (s=="pts" or s.startswith("pts_"))]
pct_col1 = list(DEdf_filt[pct_cols[0]])
pct_col2 = list(DEdf_filt[pct_cols[1]])
pct_passed = list(map(
    lambda i: (pct_col1[i] >= min_pct or pct_col2[i] >= min_pct),
    list(DEdf_filt.index)
))
DEdf_filt = DEdf_filt.loc[pct_passed]

## Output Filtered Results
DEdf_filt.to_csv(output_path('filtered_diffexp.csv'))