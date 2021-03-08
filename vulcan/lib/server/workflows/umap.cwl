cwlVersion: v1.1
class: Workflow

inputs:
  min_nCounts:
    type: int
    default: 200
    label: 'Cell Filter: minimum number of UMIs per cell'
  max_nCounts:
    type: int
    default: 30000
    label: 'Cell Filter: maximum number of UMIs per cell'
  min_nFeatures:
    type: int
    default: 100
    label: 'Cell Filter: minimum number of genes per cell'
  max_per_mito:
    type: float
    default: 20
    label: 'Cell Filter: maximum percentage of reads per cell coming from mitochondrial genes (from 0 to 100)'
  max_per_ribo:
    type: float
    default: 100
    label: 'Cell Filter: maximum percentage of reads per cell coming from ribosomal genes (from 0 to 100)'
  regress_counts:
    type: boolean
    default: true
    label: 'regress by number of counts per cell?'
  regress_genes:
    type: boolean
    default: false
    label: 'regress by number of genes per cell?'
  regress_pct_mito:
    type: boolean
    default: true
    label: 'regress by percent of reads from mitochondrial genes?'
  regress_pct_ribo:
    type: boolean
    default: false
    label: 'regress by percent of reads from ribosomal genes?'
  max_pc:
    type: int
    default: 15
    label: 'Maximum number of PCs'

outputs:
  the_data:
    type: File
    outputSource: doUmapStuff/umap

steps:
  queryMagma:
    run: scripts/fake_query.cwl
    label: 'Fetch pool record names'
    in:
      b: includeParamA
      c: includeParamB
      d: includeParamC
      e: maxPcs
    out: [names]
  pickPools:
    run: ui-queries/select-autocomplete.cwl
    label: 'Select pool records'
    in:
      a: queryMagma/names
    out: [names]
  magma_query_paths:
    run: scripts/magma_query_paths.cwl
    label: 'Retrieve path to raw counts files'
    in:
      record_ids: pickPools/names
    out: [h5_locations]
  merge_anndata_from_raw_h5:
    run: scripts/merge_anndata_from_raw_h5.cwl
    label: 'Read into scanpy and merge all records'
    in:
      h5_locations: magma_query_paths/h5_locations
    out: [merged_anndata.h5ad]
  subset_normalize_and_select_features:
    run: scripts/subset_normalize_and_select_features.cwl
    label: 'Subset cells (and then normalize)'
    in:
      merged_anndata.h5ad: merge_anndata_from_raw_h5/merged_anndata.h5ad
      min_nCounts: min_nCounts
      max_nCounts: max_nCounts
      min_nFeatures: min_nFeatures
      max_per_mito: max_per_mito
      max_per_ribo: max_per_ribo
    out: [normed_anndata.h5ad]
  regress_and_pca:
    run: scripts/regress_and_pca.cwl
    label: 'Regress parameters (and then calculate PCA)'
    in:
      normed_anndata.h5ad: subset_normalize_and_select_features/normed_anndata.h5ad
      regress_counts: regress_counts
      regress_genes: regress_genes
      regress_pct_mito: regress_pct_mito
      regress_pct_ribo: regress_pct_ribo
    out: [pca_anndata.h5ad]
  calc_umap:
    run: scripts/calc_umap.cwl
    label: 'Calculate UMAP (based on PCA)'
    in:
      pca_anndata.h5ad: regress_and_pca/pca_anndata.h5ad
      max_pc: max_pc
    out: [umap_anndata.h5ad]
  downloadRawData:
    run: ui-outputs/link.cwl
    in:
      a: calc_umap/umap_anndata.h5ad
    out: []
    label: 'Download Raw Data as h5ad'