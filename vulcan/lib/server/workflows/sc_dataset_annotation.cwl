cwlVersion: v1.1
class: Workflow

inputs:
  1_Object_Selection__dataset_name:
    type: string
    default: ''
    label: 'name of record to target'
    doc: 'Provide the identifier ("name" attribute) of the sc_seq_dataset record you wish to explore here. Note that we definitely plan to improve this particular selection method in the future!'


outputs:
  annots.csv:
    type: File
    format: csv
    outputSource: convert_map_to_file_outputs/annots.csv
  annots.xlsx:
    type: File
    format: xlsx
    outputSource: convert_map_to_file_outputs/annots.xlsx

steps:
  get_dataset_and_summarize:
    run: scripts/get_dataset_and_summarize.cwl
    label: 'Retrieve dataset & Determine metadata options'
    in:
      dataset_name: 1_Object_Selection__dataset_name
    out: [plotting_options, scdata, discrete_metadata_summary, all_opts, continuous_opts, discrete_opts, reduction_opts]
  prep_clustering_selection_options:
    run: scripts/further_summarize_metadata.cwl
    label: 'Further metadata summary'
    in:
      plotting_options: get_dataset_and_summarize/plotting_options
      discrete_metadata_summary: get_dataset_and_summarize/discrete_metadata_summary
    out: [opts, rec]
  select_clustering:
    run: ui-queries/select-autocomplete.cwl
    label: 'Pick Clustering Metadata to Annotate'
    doc: 'Select the name of the metadata column holding the clustering that you wish to annotate. Because there is no required way that this metadatas must be named, options given in this dropdown will include all non-numeric cell metadata columns of the dataset.'
    in:
      a: prep_clustering_selection_options/opts
      recommendation: prep_clustering_selection_options/rec
    out: [cluster_meta]
  prep_annotation_input:
    run: scripts/make_blank_annots.cwl
    label: 'Prep input for Annotation UI'
    in:
      cluster_meta: select_clustering/cluster_meta
      discrete_metadata_summary: get_dataset_and_summarize/discrete_metadata_summary
    out: [blank_annots.json]
  cluster_annotation:
    run: ui-queries/annotation-editor.cwl
    label: 'Annotate Clusters'
    doc: "Click the 'Edit Annotations' button, then set annotations by adding columns in the Excel-like interface. Right click in the interface for an interactions menu where you can add/remove columns. Start columns names with 'annot' to mark them as such. All other columns will be treated as for human use only, and ignored by our systems."
    in:
      a: prep_annotation_input/blank_annots.json
    out: [formulaic_data, calculated_data]
  convert_map_to_file_outputs:
    run: scripts/convert_annots_json_to_files.cwl
    label: 'Generate Output Files'
    in:
      annots.json: cluster_annotation/calculated_data
    out: [annots.csv, annots.xlsx]
  download_annots_csv:
    run: ui-outputs/link.cwl
    in:
      a: convert_map_to_file_outputs/annots.csv
    out: []
    label: 'Download Annotation Map as csv'
  download_annots_xlsx:
    run: ui-outputs/link.cwl
    in:
      a: convert_map_to_file_outputs/annots.xlsx
    out: []
    label: 'Download Annotation Map as xlsx'
  determine_new_plotting_options:
    run: scripts/add_annots_to_plotting_options.cwl
    label: 'Determine updated plotting options'
    in:
      annots.csv: convert_map_to_file_outputs/annots.csv
      plotting_options: get_dataset_and_summarize/plotting_options
      discrete_metadata_summary: get_dataset_and_summarize/discrete_metadata_summary
      discrete_opts: get_dataset_and_summarize/discrete_opts
      all_opts: get_dataset_and_summarize/all_opts
    out: [discrete_metadata_summary, discrete_opts, all_opts, plotting_options]
  plot_setup:
    run: ui-queries/any-dittoseq.cwl
    label: 'Set plot options'
    doc: "You will find your annotation columns among dropdowns that accept discrete data columns, either visible directly or visible within the 'Cell_Metadata' options set.\n\nOptions here determine both what type of plot to make, and how to set that plot up. For addtional details on individual inputs, see the 'Inputs of the Plot Configuration Interface' section of Vulcan's 'Help'-page documentation OR dittoSeq's own documentation. Any inputs without an exact label match within that Vulcan 'Help'-page table will map directly to some dittoSeq input for the function with the same name as your chosen plot type. You can find dittoSeq's documentation from either within R itself, run `?dittoSeq::<visualization-name>`, or from the RDRR.io page that pops up when you google the package."
    in:
      data_frame: determine_new_plotting_options/discrete_metadata_summary
      continuous_cols: get_dataset_and_summarize/continuous_opts
      discrete_cols: determine_new_plotting_options/discrete_opts
      all_cols: determine_new_plotting_options/all_opts
      reduction_opts: get_dataset_and_summarize/reduction_opts
    out: [plot_setup]
  make_plot:
    run: scripts/make_dittoSeq_plot_after_pulling_annots.cwl
    label: 'Create Plot'
    in:
      annots.csv: convert_map_to_file_outputs/annots.csv
      plot_setup: plot_setup/plot_setup
      scdata: get_dataset_and_summarize/scdata
      plotting_options: determine_new_plotting_options/plotting_options
    out: [plot.out, legend.png, plot.png, plot.Rds]
  download_plot:
    run: ui-outputs/link.cwl
    in:
      a: make_plot/plot.Rds
    out: []
    label: 'Download plot object'
  show_plot:
    run: ui-outputs/plot.cwl
    label: 'Display Plot'
    in:
      a: make_plot/plot.out
    out: []