cwlVersion: v1.1
class: Workflow

inputs:
  1_Cell_Filtering__min_nCounts:
    type: int
    default: 200
    label: 'minimum UMIs'
    doc: 'Minimum number of unique transcript reads per cell. Cells with fewer will be discarded. Recommended range: 200 (lower quality data, macrophage retention) to 1000 (very high quality data). Set to 0 to skip.'
  1_Cell_Filtering__max_nCounts:
    type: int
    default: 30000
    label: 'maximum UMIs'
    doc: 'Maximum number of unique transcript reads per cell. Cells with more will be discarded. Can serve as a poor-mans doublet filter. Recommended range: at least 20000. Set absurdly high (e.g. 1000000) to effectively turn this off.'
  1_Cell_Filtering__min_nFeatures:
    type: int
    default: 100
    label: 'minimum genes'
    doc: 'Minimum number of genes captured per cell. Cells with fewer will be discarded. Reasonably should be be at least 100 to 500+ for highest quality data. Set to 0 to skip.'
  1_Cell_Filtering__max_per_mito:
    type: float
    default: 20
    label: 'maximum percent mitochondrial'
    doc: 'Maximum percentage of reads per cell coming from mitochondrial genes (from 0 to 100). Cells with higher than the given percentage will be discarded. High mitochondrial content is thought to be a sign of dead/dying/low quality cells. Recommended: 5-20 depending on the data.'
  1_Cell_Filtering__max_per_ribo:
    type: float
    default: 100
    label: 'maximum percent ribosomal'
    doc: 'Maximum percentage of reads per cell coming from ribosomal genes (from 0 to 100). Cells with higher than the given percentage will be discarded. High ribosomal content essentilaly just means less "meaningful" content but is not thought to necessarily mark lower quality cells. Recommended: 100 = off or 50+ when used.'
  2_Regress_by__regress_counts:
    type: boolean
    default: true
    label: 'regress Counts?'
    doc: 'Controls whether to regress data that correlates with cells UMI counts. Regression of data confounders for PCA/UMAP/clustering calculations can improve results of these steps. NOTE: You should only regress by counts or genes, but not both'
  2_Regress_by__regress_genes:
    type: boolean
    default: false
    label: 'regress Genes?'
    doc: 'Controls whether to regress data that correlates with cells genes counts. Regression of data confounders for PCA/UMAP/clustering calculations can improve results of these steps. NOTE: You should only regress by counts or genes, but not both'
  2_Regress_by__regress_pct_mito:
    type: boolean
    default: true
    label: 'regress percent.mitochondrial?'
    doc: 'Controls whether to regress data that correlates with cells percentage of mitochondrial reads. Regression of data confounders for PCA/UMAP/clustering calculations can improve results of these steps.'
  2_Regress_by__regress_pct_ribo:
    type: boolean
    default: false
    label: 'regress percent.ribosomal?'
    doc: 'Controls whether to regress data that correlates with cells percentage of ribosomal reads. Regression of data confounders for PCA/UMAP/clustering calculations can improve results of these steps.'
  2_Regress_by__regress_tube_id:
    type: boolean
    default: false
    label: 'regress on tube IDs?'
    doc: 'Controls whether to regress data that correlates with tube IDs. Regression by this data for PCA/UMAP/clustering can be an effective form of batch correction for these steps. (We do plan to add additional batch correction options in the future!)'
  3_For_UMAP_and_Clustering__max_pc:
    type: int
    default: 15
    label: 'Number of PCs to use'
    doc: 'Principal components, from 1 to this number will be carried forward into UMAP and clustering calculations. Commonly, some number 15 or fewer is ideal. Additional tooling is planned to power tuning this parameter.'
    label: 'Maximum number of PCs'
  3_For_UMAP_and_Clustering__n_neighbors:
    type: int
    default: 10
    label: 'Number of neighbors?'
    doc: 'The size of the local neighborhood used in manifold approximation. Larger values result in a more global view, smaller values a more local view.'
  5_Cluster_Calculation__leiden_resolution:
    type: float
    default: 1
    label: 'Resolution of Leiden clustering'
    doc: 'Higher resolution yields more clusters'
  5_Cluster_Calculation__leiden_use_weights:
    type: boolean
    default: true
    label: 'use weights in leiden?'
    doc: 'If true, edge weights in the nearest neighbors graph will be used to calculate clusters.'
  4_UMAP_Calculation__umap_spread:
    type: float
    default: 1.0
    label: 'UMAP point spread?'
    doc: 'The scale of the embedded points - larger values move points further apart.'
  4_UMAP_Calculation__umap_min_dist:
    type: float
    default: 0.5
    label: 'Minimum distance?'
    doc: 'The minimum distance between cluster points - larger values produce less clumping.'
  4_UMAP_Calculation__umap_num_iters:
    type: int
    label: 'Number of iterations?'
    default: 0
    doc: 'The number of iterations for optimization - by default (0) either 200 for small datasets or 500 for large ones.'
  6_Cluster_Differential_Expression__ignore_prefixes:
    type: string
    default: 'MT-,RPL,RPS'
    label: 'Gene prefixes to ignore'
    doc: 'A set of strings, separated by commas, for which gene symbols starting with these strings should NOT be shown in the umap overlay. Not case-sensitive, so the same strings will work for both human and mouse ribo & mito genes. Note: this does not affect the full differential expression table that can be downloaded.'
  6_Cluster_Differential_Expression__dge_method:
    type: string
    default: 'wilcoxon'
    label: 'testing method'
    doc: 'A string indicating what scanpy method option to use for calculating differential expression. Options are: "logreg", "t-test", "wilcoxon", "t-test_overestim_var". See documentation for "scanpy.tl.rank_genes_groups" for further details.'

outputs:
  the_data:
    type: File
    outputSource: Differential_Expression__between_clusters/umap_workflow_anndata.h5ad

steps:
  projectData:
    run: scripts/umap-projects.cwl
    label: 'Fetch project settings'
    in: []
    out: [project_data]
  queryMagma:
    run: scripts/retrieve_selection_options.cwl
    label: 'Fetch selection options'
    in:
      project_data: projectData/project_data
    out: [selection_options]
  selectOnFeatures:
    run: ui-queries/multiple-multiselect-string-all.cwl
    label: 'Record Selection'
    doc: 'Selections here pick the subset of tube records to process and analyze. Select the values of the given features that you would like to target. The union of single-cell tube records that meet these criteria will be presented for confirmation, in the next step, based on the union of ALL feature selections here. Selections must be made for each in order to proceed, but if you want to just select tube records directly, pick the `All` option for all dropdowns here.'
    in:
      a: queryMagma/selection_options
    out: [selected_options]
  parse_record_selections:
    run: scripts/parse_record_selections.cwl
    label: 'Interpret record selection inputs.'
    in:
      project_data: projectData/project_data
      selected_options: selectOnFeatures/selected_options
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
      project_data: projectData/project_data
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
      min_nCounts: 1_Cell_Filtering__min_nCounts
      max_nCounts: 1_Cell_Filtering__max_nCounts
      min_nFeatures: 1_Cell_Filtering__min_nFeatures
      max_per_mito: 1_Cell_Filtering__max_per_mito
      max_per_ribo: 1_Cell_Filtering__max_per_ribo
    out: [normed_anndata.h5ad]
  regress_and_pca:
    run: scripts/regress_and_pca.cwl
    label: 'Regress params and run PCA'
    in:
      normed_anndata.h5ad: subset_normalize_and_select_features/normed_anndata.h5ad
      regress_counts: 2_Regress_by__regress_counts
      regress_genes: 2_Regress_by__regress_genes
      regress_pct_mito: 2_Regress_by__regress_pct_mito
      regress_pct_ribo: 2_Regress_by__regress_pct_ribo
      regress_tube_id: 2_Regress_by__regress_tube_id
    out: [pca_anndata.h5ad]
  neighbors:
    run: scripts/neighbors.cwl
    label: 'Calculate nearest neighbors (based on PCA)'
    in:
      pca_anndata.h5ad: regress_and_pca/pca_anndata.h5ad
      max_pc: 3_For_UMAP_and_Clustering__max_pc
      n_neighbors: 3_For_UMAP_and_Clustering__n_neighbors
    out: [nn_anndata.h5ad]
  calc_umap:
    run: scripts/calc_umap.cwl
    label: 'Calculate UMAP'
    in:
      spread: 4_UMAP_Calculation__umap_spread
      min_dist: 4_UMAP_Calculation__umap_min_dist
      num_iters: 4_UMAP_Calculation__umap_num_iters
      nn_anndata.h5ad: neighbors/nn_anndata.h5ad
    out: [umap_anndata.h5ad]
  calc_leiden:
    run: scripts/calc_leiden.cwl
    label: 'Calculate Leiden clustering'
    in:
      project_data: projectData/project_data
      umap_anndata.h5ad: calc_umap/umap_anndata.h5ad
      leiden_resolution: 5_Cluster_Calculation__leiden_resolution
      use_weights: 5_Cluster_Calculation__leiden_use_weights
    out: [leiden.json,leiden_anndata.h5ad,blank_annots.json]
  cluster_annotation:
    run: ui-queries/multiple-string.cwl
    label: 'Manually Annotate Clusters'
    in:
      a: calc_leiden/blank_annots.json
    out: [annots.json]
  determine_color_by_options:
    run: scripts/determine_color_by_options.cwl
    label: 'Determine coloring options'
    in:
      project_data: projectData/project_data
      leiden_anndata.h5ad: calc_leiden/leiden_anndata.h5ad
    out: [color_options, hide]
  select_color_by_option:
    run: ui-queries/scatter-plotly.cwl
    label: 'Color Options'
    in:
      data_options: determine_color_by_options/color_options
      hide: determine_color_by_options/hide
    out: [plot_setup]
  plot_umap:
    run: scripts/plot_umap.cwl
    label: 'Create UMAP plot'
    in:
      project_data: projectData/project_data
      umap_anndata.h5ad: calc_umap/umap_anndata.h5ad
      leiden.json: calc_leiden/leiden.json
      annots.json: cluster_annotation/annots.json
      color_by: select_color_by_option/plot_setup
      top10.json: Differential_Expression__between_clusters/top10.json
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
      a: Differential_Expression__between_clusters/umap_workflow_anndata.h5ad
    out: []
    label: 'Download data as h5ad'
  Differential_Expression__between_clusters:
    run: scripts/DE_btwn_clusters.cwl
    label: 'Diff. Exp.: Cluster Markers'
    in:
      leiden_anndata.h5ad: calc_leiden/leiden_anndata.h5ad
      ignore_prefixes: 6_Cluster_Differential_Expression__ignore_prefixes
      dge_method: 6_Cluster_Differential_Expression__dge_method
    out: [umap_workflow_anndata.h5ad, diffexp.csv,top10.json]
  downloadDEData:
    run: ui-outputs/link.cwl
    in:
      a: Differential_Expression__between_clusters/diffexp.csv
    out: []
    label: 'Download cluster DiffExp as csv'

