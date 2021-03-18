cwlVersion: v1.1
class: Workflow

inputs:
  min_nCounts:
    type: int
    default: 200
    label: 'minimum UMIs per cell'
    group: 'Cell Filter'
    doc: 'Minimum number of unique transcript reads per cell. Cells with fewer will be discarded. Recommended range: 200 (lower quality data, macrophage retention) to 1000 (very high quality data). Set to 0 to skip.'
  max_nCounts:
    type: int
    default: 30000
    label: 'maximum UMIs per cell'
    group: 'Cell Filter'
    doc: 'Maximum number of unique transcript reads per cell. Cells with more will be discarded. Can serve as a poor-mans doublet filter. Recommended range: at least 20000. Set absurdly high (e.g. 1000000) to effectively turn this off.'
  min_nFeatures:
    type: int
    default: 100
    label: 'minimum genes per cell'
    group: 'Cell Filter'
    doc: 'Minimum number of genes captured per cell. Cells with fewer will be discarded. Reasonably should be be at least 100 to 500+ for highest quality data. Set to 0 to skip.'
  max_per_mito:
    type: float
    default: 20
    label: 'maximum percent mitochondrial'
    group: 'Cell Filter'
    doc: 'Maximum percentage of reads per cell coming from mitochondrial genes (from 0 to 100). Cells with higher than the given percentage will be discarded. High mitochondrial content is thought to be a sign of dead/dying/low quality cells. Recommended: 5-20 depending on the data.'
  max_per_ribo:
    type: float
    default: 100
    label: 'maximum percent ribosomal'
    group: 'Cell Filter'
    doc: 'Maximum percentage of reads per cell coming from ribosomal genes (from 0 to 100). Cells with higher than the given percentage will be discarded. High ribosomal content essentilaly just means less "meaningful" content but is not thought to necessarily mark lower quality cells. Recommended: 100 = off or 50+ when used.'
  regress_counts:
    type: boolean
    default: true
    label: 'regress Counts?'
    group: 'Regress'
    doc: 'Controls whether to regress data that correlates with cells UMI counts. Regression of data confounders prior to PCA/UMAP/clustering can improve results of these steps. NOTE: You should only regress by counts or genes, but not both'
  regress_genes:
    type: boolean
    default: false
    label: 'regress Genes?'
    group: 'Regress'
    doc: 'Controls whether to regress data that correlates with cells genes counts. Regression of data confounders prior to PCA/UMAP/clustering can improve results of these steps. NOTE: You should only regress by counts or genes, but not both'
  regress_pct_mito:
    type: boolean
    default: true
    label: 'regress percent.mitochondrial?'
    group: 'Regress'
    doc: 'Controls whether to regress data that correlates with cells percentage of mitochondrial reads. Regression of data confounders prior to PCA/UMAP/clustering can improve results of these steps.'
  regress_pct_ribo:
    type: boolean
    default: false
    label: 'regress percent.ribosomal?'
    group: 'Regress'
    doc: 'Controls whether to regress data that correlates with cells percentage of ribosomal reads. Regression of data confounders prior to PCA/UMAP/clustering can improve results of these steps.'
  max_pc:
    type: int
    default: 15
    label: 'Number of PCs to use'
    group: 'UMAP Calculation'
    doc: 'Principal components, from 1 to this number will be carried forward into UMAP and clustering calculations. Commonly, some number 15 or fewer is ideal. Additional tooling is planned to power tuning this parameter.'

outputs:
  the_data:
    type: File
    outputSource: calc_umap/umap_anndata.h5ad

steps:
  queryMagma:
    run: scripts/fake_query.cwl
    label: 'Fetch pool record names'
    in:
      a: min_nCounts
      b: max_nCounts
      c: min_nFeatures
      d: max_per_mito
      e: max_per_ribo
      f: regress_counts
      g: regress_genes
      h: regress_pct_mito
      i: regress_pct_ribo
      j: max_pc
    out: [names]
  pickPools:
    run: ui-queries/select-autocomplete.cwl
    label: 'Select pool records'
    in:
      a: queryMagma/names
    out: [names]
  magma_query_paths:
    run: scripts/magma_query_paths.cwl
    label: 'Obtain raw data locations'
    in:
      record_ids: pickPools/names
    out: [h5_locations]
  merge_anndata_from_raw_h5:
    run: scripts/merge_anndata_from_raw_h5.cwl
    label: 'Import into scanpy'
    in:
      h5_locations: magma_query_paths/h5_locations
    out: [merged_anndata.h5ad]
  subset_normalize_and_select_features:
    run: scripts/subset_normalize_and_select_features.cwl
    label: 'Subset cells and normalize'
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
    label: 'Regress params and run PCA'
    in:
      normed_anndata.h5ad: subset_normalize_and_select_features/normed_anndata.h5ad
      regress_counts: regress_counts
      regress_genes: regress_genes
      regress_pct_mito: regress_pct_mito
      regress_pct_ribo: regress_pct_ribo
    out: [pca_anndata.h5ad]
  calc_umap:
    run: scripts/calc_umap.cwl
    label: 'Calculate UMAP'
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