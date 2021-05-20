cwlVersion: v1.1
class: Workflow

inputs:
  1_Processed_data__link:
    type: str
    default: ''
    label: 'Provide a link?'
    doc: 'TODO'

outputs:
  the_data:
    type: File
    outputSource: dge/processed_anndata.h5ad

steps:
  get_data:
    run: scrpts/get_data.cwl
    in:
      url: 1_Processed_data__link
    out: [ processed_data.h5ad ]



  processedData:
  run: scripts/get_data.cwl
    label: 'Fetch processed h5'
    in: []
    out: [processed_h5]
  queryMagma:
    run: scripts/retrieve_selection_options.cwl
    label: 'Fetch selection options'
    in:
      project_data: projectData/project_data
    out: [selection_options]
  selectOnFeatures:
    run: ui-queries/multiple-multiselect-string-all.cwl
    label: 'Record Selection'
    doc: 'Selections here pick the subset of tube records to process and analyze.
    Select the values of the given features that you would like to target.
    The union of single-cell tube records that meet these criteria will be presented for confirmation, in the next step,
    based on the union of ALL feature selections here. Selections must be made for each in order to proceed,
    but if you want to just select tube records directly, pick the `All` option for all dropdowns here.'
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
      project_data: projectData/project_data
      spread: 4_UMAP_Calculation__umap_spread
      min_dist: 4_UMAP_Calculation__umap_min_dist
      num_iters: 4_UMAP_Calculation__umap_num_iters
      nn_anndata.h5ad: neighbors/nn_anndata.h5ad
    out: [umap_anndata.h5ad,color_options]
  calc_leiden:
    run: scripts/calc_leiden.cwl
    label: 'Calculate Leiden clustering'
    in:
      nn_anndata.h5ad: neighbors/nn_anndata.h5ad
      leiden_resolution: 5_Cluster_Calculation__leiden_resolution
      use_weights: 5_Cluster_Calculation__leiden_use_weights
    out: [leiden.json]
  select_color_by_option:
    run: ui-queries/nested-select-autocomplete.cwl
    label: 'Color Options'
    in:
      a: calc_umap/color_options
    out: [color_by]
  plot_umap:
    run: scripts/plot_umap.cwl
    label: 'Create UMAP plot'
    in:
      project_data: projectData/project_data
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


