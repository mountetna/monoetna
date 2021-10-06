cwlVersion: v1.1
class: Workflow

inputs:
  1_DiffExp_Calculation_Inputs__test_method:
    type: string
    default: t-test
    label: 'p-value calc method'
    doc: 'One of [logreg, t-test, wilcoxon, t-test_overestim_var]. Sets the `method` input of the `scanpy.tl.rank_genes_groups()` function used for differential expression calculation.'
  2_After_Calculation_Cutoffs__min_abs_fc:
    type: float
    default: 0
    label: 'Min. Fold Change (absolute value)'
    doc: 'A number, used for subsetting the differential expression output based on fold-changes. Common values: 0.5, 1, 2.'
  2_After_Calculation_Cutoffs__max_pval:
    type: float
    default: 0.05
    label: 'Max. Adjusted P-value'
    doc: 'A number, used for subsetting the differential expression output based on p-values.'
  2_After_Calculation_Cutoffs__min_pct:
    type: float
    default: 0.05
    label: 'Min. Percent Expression'
    doc: 'A number, used for subsetting the differential expression output based on percent capture of genes within the comparison groups. Genes are retained as long as either side of the comparison has a higher percent capture.'
  2_After_Calculation_Cutoffs__pos_only:
    type: boolean
    default: false
    label: 'Upregulated only?'
    doc: 'true or false, used for subsetting the differential expression output based fold change directionality. When true, only genes upregulated in the group-1 will be retained.'

outputs:
  the_csv:
    type: File
    outputSource: DGEcalc/diffexp.csv

steps:
  mockSetup:
    run: scripts/dge_mock_data.cwl
    label: 'Fetch data'
    in: []
    out: [scdata.h5ad]
    
  extract_metadata:
    run: scripts/dge_grab_obs.cwl
    label: 'Extract metadata Options'
    in:
      scdata.h5ad: mockSetup/scdata.h5ad
    out: [metadata]
    
  setDEMethods:
    run: ui-queries/diff-exp-sc.cwl
    label: 'Set how DE should be run'
    in:
      a: extract_metadata/metadata
    out: [dge_setup]
  
  DGEcalc:
    run: scripts/dge_calc.cwl
    label: 'Calculate DGE'
    in:
      scdata.h5ad: mockSetup/scdata.h5ad
      setup: setDEMethods/dge_setup
      test_method: 1_DiffExp_Calculation_Inputs__test_method
    out: [diffexp.csv]
    
  downloadDEData:
    run: ui-outputs/link.cwl
    in:
      a: DGEcalc/diffexp.csv
    out: []
    label: 'Download Full DiffExp Results as csv'
  
  DGEfilter:
    run: scripts/dge_filter.cwl
    label: 'Filter the DGE Output'
    in:
      full_diffexp.csv: DGEcalc/diffexp.csv
      min_abs_fc: 2_After_Calculation_Cutoffs__min_abs_fc
      max_pval: 2_After_Calculation_Cutoffs__max_pval
      min_pct: 2_After_Calculation_Cutoffs__min_pct
      pos_only: 2_After_Calculation_Cutoffs__pos_only
    out: [filtered_diffexp.csv]
  
  downloadFilteredDEData:
    run: ui-outputs/link.cwl
    in:
      a: DGEfilter/filtered_diffexp.csv
    out: []
    label: 'Download Filtered DiffExp Results as csv'
