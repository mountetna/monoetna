cwlVersion: v1.1
class: Workflow

inputs:
  CellFilter__min_nCounts:
    type: int
    default: 200
    label: 'minimum UMIs per cell'
    doc: 'Minimum number of unique transcript reads per cell. Cells with fewer will be discarded. Recommended range: 200 (lower quality data, macrophage retention) to 1000 (very high quality data). Set to 0 to skip.'
  CellFilter__max_nCounts:
    type: int
    default: 30000
    label: 'maximum UMIs per cell'
    group: 'Cell Filter'
    doc: 'Maximum number of unique transcript reads per cell. Cells with more will be discarded. Can serve as a poor-mans doublet filter. Recommended range: at least 20000. Set absurdly high (e.g. 1000000) to effectively turn this off.'
  CellFilter__min_nFeatures:
    type: int
    default: 100
    label: 'minimum genes per cell'
    doc: 'Minimum number of genes captured per cell. Cells with fewer will be discarded. Reasonably should be be at least 100 to 500+ for highest quality data. Set to 0 to skip.'
  CellFilter__max_per_mito:
    type: float
    default: 20
    label: 'maximum percent mitochondrial'
    doc: 'Maximum percentage of reads per cell coming from mitochondrial genes (from 0 to 100). Cells with higher than the given percentage will be discarded. High mitochondrial content is thought to be a sign of dead/dying/low quality cells. Recommended: 5-20 depending on the data.'
  CellFilter__max_per_ribo:
    type: float
    default: 100
    label: 'maximum percent ribosomal'
    doc: 'Maximum percentage of reads per cell coming from ribosomal genes (from 0 to 100). Cells with higher than the given percentage will be discarded. High ribosomal content essentilaly just means less "meaningful" content but is not thought to necessarily mark lower quality cells. Recommended: 100 = off or 50+ when used.'
  Regress__regress_counts:
    type: boolean
    default: true
    label: 'regress Counts?'
    doc: 'Controls whether to regress data that correlates with cells UMI counts. Regression of data confounders prior to PCA/UMAP/clustering can improve results of these steps. NOTE: You should only regress by counts or genes, but not both'
  Regress__regress_genes:
    type: boolean
    default: false
    label: 'regress Genes?'
    doc: 'Controls whether to regress data that correlates with cells genes counts. Regression of data confounders prior to PCA/UMAP/clustering can improve results of these steps. NOTE: You should only regress by counts or genes, but not both'
  Regress__regress_pct_mito:
    type: boolean
    default: true
    label: 'regress percent.mitochondrial?'
    doc: 'Controls whether to regress data that correlates with cells percentage of mitochondrial reads. Regression of data confounders prior to PCA/UMAP/clustering can improve results of these steps.'
  Regress__regress_pct_ribo:
    type: boolean
    default: false
    label: 'regress percent.ribosomal?'
    doc: 'Controls whether to regress data that correlates with cells percentage of ribosomal reads. Regression of data confounders prior to PCA/UMAP/clustering can improve results of these steps.'
  UMAP_Calculation__max_pc:
    type: int
    default: 15
    label: 'Number of PCs to use'
    doc: 'Principal components, from 1 to this number will be carried forward into UMAP and clustering calculations. Commonly, some number 15 or fewer is ideal. Additional tooling is planned to power tuning this parameter.'

outputs:
  the_data:
    type: File
    outputSource: calc_umap/umap_anndata.h5ad

steps:
  queryMagma:
    run: scripts/fake_query.cwl
    label: 'Fetch user selection options'
    in:
      a: CellFilter__min_nCounts
      b: CellFilter__max_nCounts
      c: CellFilter__min_nFeatures
      d: CellFilter__max_per_mito
      e: CellFilter__max_per_ribo
      f: Regress__regress_counts
      g: Regress__regress_genes
      h: Regress__regress_pct_mito
      i: Regress__regress_pct_ribo
      j: UMAP_Calculation__max_pc
    out: [experiments, tissues, records, pools]
  pickExperiments:
    run: ui-queries/multiselect-string.cwl
    label: 'Select experiments'
    in:
      a: queryMagma/experiments
    out: [names]
  pickTissues:
    run: ui-queries/multiselect-string.cwl
    label: 'Select tissue types'
    in:
      a: queryMagma/tissues
    out: [names]
  pickPools:
    run: ui-queries/multiselect-string.cwl
    label: 'Select pool records'
    in:
      a: queryMagma/pools
    out: [names]
  pickTubes:
    run: ui-queries/multiselect-string.cwl
    label: 'Select individual tube records'
    in:
      a: queryMagma/records
    out: [names]
  parseRecordSelections:
    run: scripts/parse_record_selections.cwl
    label: 'Interpret record selection inputs.'
    in:
      experiments:  pickExperiments/names
      tissues: pickTissues/names
      pools: pickPools/names
      tubes: pickTubes/names
    out: [tube_recs]
  verifyRecordNames:
    run: ui-queries/checkboxes.cwl
    label: 'Confirm record names'
    in:
      a: parseRecordSelections/tube_recs
    out: [names]
  magma_query_paths:
    run: scripts/magma_query_paths.cwl
    label: 'Query paths to raw counts files'
    in:
      record_ids: verifyRecordNames/names
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
      min_nCounts: CellFilter__min_nCounts
      max_nCounts: CellFilter__max_nCounts
      min_nFeatures: CellFilter__min_nFeatures
      max_per_mito: CellFilter__max_per_mito
      max_per_ribo: CellFilter__max_per_ribo
    out: [normed_anndata.h5ad]
  regress_and_pca:
    run: scripts/regress_and_pca.cwl
    label: 'Regress params and run PCA'
    in:
      normed_anndata.h5ad: subset_normalize_and_select_features/normed_anndata.h5ad
      regress_counts: Regress__regress_counts
      regress_genes: Regress__regress_genes
      regress_pct_mito: Regress__regress_pct_mito
      regress_pct_ribo: Regress__regress_pct_ribo
    out: [pca_anndata.h5ad]
  calc_umap:
    run: scripts/calc_umap.cwl
    label: 'Calculate UMAP'
    in:
      pca_anndata.h5ad: regress_and_pca/pca_anndata.h5ad
      max_pc: UMAP_Calculation__max_pc
    out: [umap_anndata.h5ad]
  downloadRawData:
    run: ui-outputs/link.cwl
    in:
      a: calc_umap/umap_anndata.h5ad
    out: []
    label: 'Download Raw Data as h5ad'