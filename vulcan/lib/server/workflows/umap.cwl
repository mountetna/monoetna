cwlVersion: v1.1
class: Workflow

inputs:
  Cell_Filtering__min_nCounts:
    type: int
    default: 200
    label: 'minimum UMIs per cell'
    doc: 'Minimum number of unique transcript reads per cell. Cells with fewer will be discarded. Recommended range: 200 (lower quality data, macrophage retention) to 1000 (very high quality data). Set to 0 to skip.'
  Cell_Filtering__max_nCounts:
    type: int
    default: 30000
    label: 'maximum UMIs per cell'
    doc: 'Maximum number of unique transcript reads per cell. Cells with more will be discarded. Can serve as a poor-mans doublet filter. Recommended range: at least 20000. Set absurdly high (e.g. 1000000) to effectively turn this off.'
  Cell_Filtering__min_nFeatures:
    type: int
    default: 100
    label: 'minimum genes per cell'
    doc: 'Minimum number of genes captured per cell. Cells with fewer will be discarded. Reasonably should be be at least 100 to 500+ for highest quality data. Set to 0 to skip.'
  Cell_Filtering__max_per_mito:
    type: float
    default: 20
    label: 'maximum percent mitochondrial'
    doc: 'Maximum percentage of reads per cell coming from mitochondrial genes (from 0 to 100). Cells with higher than the given percentage will be discarded. High mitochondrial content is thought to be a sign of dead/dying/low quality cells. Recommended: 5-20 depending on the data.'
  Cell_Filtering__max_per_ribo:
    type: float
    default: 100
    label: 'maximum percent ribosomal'
    doc: 'Maximum percentage of reads per cell coming from ribosomal genes (from 0 to 100). Cells with higher than the given percentage will be discarded. High ribosomal content essentilaly just means less "meaningful" content but is not thought to necessarily mark lower quality cells. Recommended: 100 = off or 50+ when used.'
  Regress__regress_counts:
    type: boolean
    default: true
    label: 'regress Counts?'
    doc: 'Controls whether to regress data that correlates with cells UMI counts. Regression of data confounders for PCA/UMAP/clustering calculations can improve results of these steps. NOTE: You should only regress by counts or genes, but not both'
  Regress__regress_genes:
    type: boolean
    default: false
    label: 'regress Genes?'
    doc: 'Controls whether to regress data that correlates with cells genes counts. Regression of data confounders for PCA/UMAP/clustering calculations can improve results of these steps. NOTE: You should only regress by counts or genes, but not both'
  Regress__regress_pct_mito:
    type: boolean
    default: true
    label: 'regress percent.mitochondrial?'
    doc: 'Controls whether to regress data that correlates with cells percentage of mitochondrial reads. Regression of data confounders for PCA/UMAP/clustering calculations can improve results of these steps.'
  Regress__regress_pct_ribo:
    type: boolean
    default: false
    label: 'regress percent.ribosomal?'
    doc: 'Controls whether to regress data that correlates with cells percentage of ribosomal reads. Regression of data confounders for PCA/UMAP/clustering calculations can improve results of these steps.'
  Regress__regress_tube_id:
    type: boolean
    default: false
    label: 'regress on tube IDs?'
    doc: 'Controls whether to regress data that correlates with tube IDs. Regression by this data for PCA/UMAP/clustering can be an effective form of batch correction for these steps. (We do plan to add additional batch correction options in the future!)'
  UMAP_Calculation__max_pc:
    type: int
    default: 15
    label: 'Number of PCs to use'
    doc: 'Principal components, from 1 to this number will be carried forward into UMAP and clustering calculations. Commonly, some number 15 or fewer is ideal. Additional tooling is planned to power tuning this parameter.'
    label: 'Maximum number of PCs'
  UMAP_Calculation__n_neighbors:
    type: int
    default: 10
    label: 'Number of neighbors?'
    doc: 'The size of the local neighborhood used in manifold approximation. Larger values result in a more global view, smaller values a more local view.'
  UMAP_Calculation__leiden_resolution:
    type: float
    default: 1
    label: 'Resolution of Leiden clustering'
    doc: 'Higher resolution yields more clusters'
  UMAP_Calculation__leiden_use_weights:
    type: boolean
    default: true
    label: 'use weights in leiden?'
    doc: 'If true, edge weights in the nearest neighbors graph will be used to calculate clusters.'
  UMAP_Calculation__umap_spread:
    type: float
    default: 1.0
    label: 'UMAP point spread?'
    doc: 'The scale of the embedded points - larger values move points further apart.'
  UMAP_Calculation__umap_min_dist:
    type: float
    default: 0.5
    label: 'Minimum distance?'
    doc: 'The minimum distance between cluster points - larger values produce less clumping.'
  UMAP_Calculation__umap_num_iters:
    type: int
    label: 'Number of iterations?'
    default: 0
    doc: 'The number of iterations for optimization - by default (0) either 200 for small datasets or 500 for large ones.'

outputs:
  the_data:
    type: File
    outputSource: calc_umap/umap_anndata.h5ad

steps:
  queryMagma:
    run: scripts/retrieve_selection_options.cwl
    label: 'Fetch selection options'
    in:
      a: Cell_Filtering__min_nCounts
      b: Cell_Filtering__max_nCounts
      c: Cell_Filtering__min_nFeatures
      d: Cell_Filtering__max_per_mito
      e: Cell_Filtering__max_per_ribo
      f: Regress__regress_counts
      g: Regress__regress_genes
      h: Regress__regress_pct_mito
      i: Regress__regress_pct_ribo
      j: Regress__regress_tube_id
      k: UMAP_Calculation__max_pc
      l: UMAP_Calculation__leiden_resolution
    out: [experiments, tissues, color_options]
  Select_Records__pickExperiments:
    run: ui-queries/multiselect-string.cwl
    label: 'Select Experiments'
    doc: 'Subset  of experiments to use, based on their `alias`. These selections get combined with Tissue selections with AND logic. If you want to just select tube records directly, pick No Selections for all dropdowns here.'
    in:
      a: queryMagma/experiments
    out: [options]
  Select_Records__pickTissues:
    run: ui-queries/multiselect-string.cwl
    label: 'Select Tissues'
    doc: 'Subset of biospecimen_types to use. These selections get combined with Experiment selections with AND logic. If you want to just select tube records directly, pick No Selections for all dropdowns here.'
    in:
      a: queryMagma/tissues
    out: [options]
  parse_record_selections:
    run: scripts/parse_record_selections.cwl
    label: 'Interpret record selection inputs.'
    in:
      experiments: Select_Records__pickExperiments/options
      tissues: Select_Records__pickTissues/options
    out: [tube_recs]
  verifyRecordNames:
    run: ui-queries/checkboxes.cwl
    label: 'Confirm record names'
    in:
      a: parse_record_selections/tube_recs
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
      min_nCounts: Cell_Filtering__min_nCounts
      max_nCounts: Cell_Filtering__max_nCounts
      min_nFeatures: Cell_Filtering__min_nFeatures
      max_per_mito: Cell_Filtering__max_per_mito
      max_per_ribo: Cell_Filtering__max_per_ribo
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
      regress_tube_id: Regress__regress_tube_id
    out: [pca_anndata.h5ad]
  neighbors:
    run: scripts/neighbors.cwl
    label: 'Calculate nearest neighbors (based on PCA)'
    in:
      pca_anndata.h5ad: regress_and_pca/pca_anndata.h5ad
      max_pc: UMAP_Calculation__max_pc
      n_neighbors: UMAP_Calculation__n_neighbors
    out: [nn_anndata.h5ad]
  calc_umap:
    run: scripts/calc_umap.cwl
    label: 'Calculate UMAP'
    in:
      spread: UMAP_Calculation__umap_spread
      min_dist: UMAP_Calculation__umap_min_dist
      num_iters: UMAP_Calculation__umap_num_iters
      nn_anndata.h5ad: neighbors/nn_anndata.h5ad
    out: [umap_anndata.h5ad]
  calc_leiden:
    run: scripts/calc_leiden.cwl
    label: 'Calculate Leiden clustering'
    in:
      nn_anndata.h5ad: neighbors/nn_anndata.h5ad
      leiden_resolution: UMAP_Calculation__leiden_resolution
    out: [leiden.json]
  select_color_by_option:
    run: ui-queries/nested-select-autocomplete.cwl
    label: 'Color Options'
    in:
      a: queryMagma/color_options
      b: calc_leiden/leiden.json
    out: [color_by]
  plot_umap:
    run: scripts/plot_umap.cwl
    label: 'Create UMAP plot'
    in:
      umap_anndata.h5ad: calc_umap/umap_anndata.h5ad
      leiden.json: calc_leiden/leiden.json
      color_by: select_color_by_option/color_by
    out: [umap.plotly.json]
  show_umap_plot:
    run: ui-outputs/plotly.cwl
    in:
      a: plot_umap/umap.plotly.json
    out: []
    label: 'Display UMAP'
  downloadRawData:
    run: ui-outputs/link.cwl
    in:
      a: calc_umap/umap_anndata.h5ad
    out: []
    label: 'Download data as h5ad'

