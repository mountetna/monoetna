cwlVersion: v1.1
class: Workflow

inputs:
  min_nCounts:
    type: int
    default: 200
    label: 'minimum counts per cell'
  max_nCounts:
    type: int
    default: 30000
    label: 'maximum counts per cell'
  min_nFeatures:
    type: int
    default: 100
    label: 'minimum genes/features per cell'
  max_per_mito:
    type: int
    default: 20
    label: 'maximum percent of mitochondrial reads per cell'
  max_per_ribo:
    type: int
    default: 50
    label: 'maximum percent of ribosomal reads per cell'
  regress_nCounts:
    type: boolean
    default: false
    label: nCounts
  regress_nFeatures:
    type: boolean
    default: false
    label: nFeatures
  regress_per_mito:
    type: boolean
    default: false
    label: '% mito'
  regress_per_ribo:
    type: boolean
    default: false
    label: '% ribo'
  max_pc:
    type: int
    default: 15
    label: 'Maximum number of Principal Components'

outputs:
  subset_data:
    type: File
    outputSource: subset_normalize_and_select_features/sc_rna_seq_normalized_data
  pc_data:
    type: File
    outputSource: regress_and_pca/pca_output_data
  umap_data:
    type: File
    outputSource: umap/umap_data
  expression_matrix:
    type: File
    outputSource: umap/expression_matrix

steps:
  query_magma:
    run: scripts/query-mvir1-sc_rna_seq_pools.cwl
    in: []
    out: [sc_rna_seq_pool_names]
  ui_pick_subset:
    run: ui-queries/ui-query-magma.cwl
    in:
      all_pool_names: query_magma/sc_rna_seq_pool_names
    out: [subset_sc_rna_pool_names]
  subset_normalize_and_select_features:
    run: scripts/subset_normalize_and_select_features.cwl
    in:
      sc_rna_seq_pools: ui_pick_subset/subset_sc_rna_pool_names
      min_nCounts: min_nCounts
      max_nCounts: max_nCounts
      min_nFeatures: min_nFeatures
      max_per_mito: max_per_mito
      max_per_ribo: max_per_ribo
    out: [sc_rna_seq_normalized_data]
  regress_and_pca:
    run: scripts/regress_and_pca.cwl
    in:
      data: subset_normalize_and_select_features/sc_rna_seq_normalized_data
      counts: regress_nCounts
      genes: regress_nFeatures
      per_mito: regress_per_mito
      per_ribo: regress_per_ribo
    out: [pca_output_data]
  umap:
    run: scripts/umap.cwl
    in:
      data: regress_and_pca/pca_output_data
      max_pc: max_pc
    out: [umap_data, expression_matrix]
  ui_plot:
    run: xy.cwl # The type of plot.
    in:
      series0: umap/umap_data # Specify the primary plot for x-y coordinates
      group0: umap/expression_matrix # Specify any secondary data to color by
      series0__type:
        default: scatter # The type of series
    out: []
