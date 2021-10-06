from archimedes.functions.dataflow import input_json, input_var, input_path, output_json
from archimedes.functions.scanpy import scanpy as sc
from archimedes.functions.utils import pandas as pd
from archimedes.functions.list import unique

# Functions for DE Extraction and Output
def DF_per_group(DEdat, meta_used, group):
    # keys = DEdat.keys
    DF = pd.DataFrame(
        dict([
            this_col,
            DEdat[this_col][group]
        # ] for this_col in keys
        ] for this_col in ['names', 'scores', 'pvals', 'pvals_adj', 'logfoldchanges']
        )
    )
    DF[meta_used] = group
    return DF
def DF_all_groups(DEdat, meta_used, group_use = None):
    if group_use==None:
        group_use = unique(scdata.obs['leiden'])
    return pd.concat( (DF_per_group(DEdat, meta_used, group)) for group in group_use )

# Read inputs
scdata = sc.read(input_path('scdata.h5ad'))
test_method = input_var('test_method')
setup = input_json('setup')
# And parse the minimal set of setup inputs
method = setup['method']
de_meta = setup['de_meta']

# Subset if there are subset inputs in 'setup'
if 'subset_meta' in setup.keys:
    subset_meta = setup['subset_meta']
    subset_use = setup['subset_use']
    scdata = scdata[scdata.obs[subset_meta].isin(subset_use)]

## Run DE
# Method 1:
if setup['method']=='btwn-all-de-groups':
    
    # Run DE
    sc.tl.rank_genes_groups(
        scdata, method=test_method,
        groupby=de_meta,
        # use_raw=True, layer=...,
        pct=True, reference='rest', corr_method='benjamini-hochberg')
    
    # Extract results and reshape
    DEdat = scdata.uns['rank_genes_groups']
    DEdf = DF_all_groups(DEdat, de_meta)

# Method 2:
if method=='btwn-sets':
    de_group_1 = setup['de_group_1']
    de_group_2 = setup['de_group_2']
    
    # Subset to cells in these groups
    scdata = scdata[ scdata.obs[ de_meta ].isin( de_group_1 + de_group_2 )]
    
    # Add metadata of exactly g1 vs g2 
    string_g1 = '_'.join(de_group_1)
    string_g2 = '_'.join(de_group_2)
    scdata.obs['DE_groups'] = list(
        {True: string_g1, False: string_g2}[is_in] for is_in in scdata.obs[ de_meta ].isin( de_group_1 )
    )
    
    # Run DE
    sc.tl.rank_genes_groups(
        scdata, method=test_method,
        groupby='DE_groups', groups = [string_g1 + string_g2],        
        # use_raw=True, layer=...,
        pct=True, reference='rest', corr_method='benjamini-hochberg')
    
    # Extract results and reshape
    DEdat = scdata.uns['rank_genes_groups']
    DEdf = DF_all_groups(DEdat, 'DE_groups')

# Method 3
# if method=='btwn-sets-multiple-groups':
    #

# OUTPUT.to_csv(output_path('diffexp.csv'))

## Filter

# OUTPUT.to_csv(output_path('filtered_diffexp.csv'))