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
  select_batch_options:
    run: ui-queries/select-autocomplete.cwl
    label: 'Pick Clustering Metadata to Annotation'
    doc: 'Selects the data to use in subsequent steps.'
    in:
      a: prep_clustering_selection_options/opts
      recommendation: prep_clustering_selection_options/rec
    out: [cluster_meta]
  prep_annotation_input:
    run: scripts/make_blank_annots.cwl
    label: 'Prep input for Annotation UI'
    in:
      cluster_meta: select_batch_options/cluster_meta
      discrete_metadata_summary: get_dataset_and_summarize/discrete_metadata_summary
    out: [blank_annots.json]
  cluster_annotation:
    run: ui-queries/annotation-editor.cwl
    label: 'Annotate Clusters'
    doc: "Set annotations by adding columns in the Excel-like interface. Right click for an interactions menu where you can add/remove columns. Start columns names with 'annot' to mark them as such. All other columns will be treated as for human use only, and ignored by our system."
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