from archimedes.functions.dataflow import input_json, input_var, input_path, output_path
from archimedes.functions.scanpy import scanpy as sc
from archimedes.functions.utils import pandas as pd
from archimedes.functions.list import unique, ensure_list

# Functions for DE Extraction and Output
def collect_DF_per_set(DEdat, meta_used, set, group_name = None):
    # keys = DEdat.keys()
    keys = ['names', 'pvals', 'pvals_adj', 'logfoldchanges']
    DF = pd.DataFrame(
        dict(
            [
                this_col,
                DEdat[this_col][set]
            ] for this_col in keys
        )
    )
    ## Add pts (percent expressed), ensuring proper ordering
    names = list(DF['names'])
    if 'pts_rest' in DEdat.keys():
        for key in ['pts', 'pts_rest']:
            DF[key] = list(DEdat[key].reindex(names)[set])
    if not 'pts_rest' in DEdat.keys():
        # Gotta grab and rename.
        grps = set.split("__VS__")
        pct_dat = DEdat['pts'].reindex(names)
        for i in [0,1]:
            DF['pts_' + grps[i]] = list(pct_dat[pct_dat.keys()[i]])
    DF[meta_used] = set
    if group_name!=None:
        DF['DE_group'] = group_name
    return DF
def collect_DF_all_sets(DEdat, meta_used, scdata, sets_use = None):
    if sets_use==None:
        sets_use = unique(scdata.obs[meta_used])
    return pd.concat( (collect_DF_per_set(DEdat, meta_used, set)) for set in sets_use )

# Function for running set-vs-set targeted DE
def DE_specific_sets(scdata, de_meta, de_g1, de_g2, test_method, group = None):
    # Subset to cells in these groups
    scdata_setsub = scdata[ scdata.obs[ de_meta ].isin( de_g1 + de_g2 )]
    
    # Add metadata of exactly g1 vs g2 
    string_g1 = '_'.join(de_g1)
    string_g2 = '_'.join(de_g2)
    string_report = string_g1 + "__VS__" + string_g2
    scdata_setsub.obs['DE_sets'] = list(
        {True: string_report, False: 'never_seen'}[is_in] for is_in in scdata_setsub.obs[ de_meta ].isin( de_g1 )
    )
    
    # Run DE
    sc.tl.rank_genes_groups(
        scdata_setsub, method=test_method,
        groupby='DE_sets', groups = [string_report], reference='never_seen',        
        # use_raw=True, layer=...,
        pts=True, corr_method='benjamini-hochberg')
    
    # Extract results and reshape
    DEdat = scdata_setsub.uns['rank_genes_groups']
    # return DEdat
    return collect_DF_per_set(DEdat, 'DE_sets', string_report, group_name=group)

## Read inputs
scdata = sc.read(input_path('scdata.h5ad'))
test_method = input_var('test_method')
setup = input_json('setup')
# And parse the minimal set of setup inputs
method = setup['method']
de_meta = setup['de_meta']

## Subset if there are subset inputs in 'setup'
if 'subset_meta' in setup.keys():
    subset_meta = setup['subset_meta']
    subset_use = ensure_list(setup['subset_use'])
    scdata = scdata[scdata.obs[subset_meta].isin(subset_use)]

## Run DE
# Method 1:
if setup['method']=='btwn-all-de-groups':
    
    # Run DE
    sc.tl.rank_genes_groups(
        scdata, method=test_method,
        groupby=de_meta, reference='rest',
        # use_raw=True, layer=...,
        pts=True, corr_method='benjamini-hochberg')
    
    # Extract results and reshape
    DEdat = scdata.uns['rank_genes_groups']
    DEdf = collect_DF_all_sets(DEdat, de_meta, scdata)

# Method 2:
if method=='btwn-sets':
    de_group_1 = ensure_list(setup['de_group_1'])
    de_group_2 = ensure_list(setup['de_group_2'])
    DEdf = DE_specific_sets(scdata, de_meta, de_group_1, de_group_2, test_method, None)

# Method 3
if method=='btwn-sets-multiple-groups':
    de_group_1 = ensure_list(setup['de_group_1'])
    de_group_2 = ensure_list(setup['de_group_2'])
    group_meta = setup['group_meta']
    group_use = ensure_list(setup['group_use'])
    
    subsets = dict(
        [
            group,
            scdata[scdata.obs[group_meta]==group]
        ] for group in group_use
    )
    DEdf = pd.concat(
            DE_specific_sets(subsets[group], de_meta, de_group_1, de_group_2, test_method, group = group)
            for group in group_use
        )

## Output unfiltered
DEdf.to_csv(output_path('full_diffexp.csv'))