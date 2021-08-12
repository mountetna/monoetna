cwlVersion: v1.1
class: Workflow

inputs:
  1_ignore__spacer:
    type: int
    default: 200
    label: "Ignore me"
    doc: "Does nothing"

outputs:
  the_plot:
    type: File
    outputSource: plot_scatter/scatter_plot.json

steps:
  get_data:
    run: scripts/mockDF.cwl
    label: 'Fetch Data'
    in: []
    out: [data_frame, data_options]
  fill_plot_options:
    run: ui-queries/scatter-plotly.cwl
    label: 'Set plot options'
    in:
      data_options: get_data/data_options
      data_frame: get_data/data_frame
    out: [data_options]
  plot_scatter:
    run: scripts/scatter_make_scatter.cwl
    label: 'Create Scatter Plot'
    in:
      plot_settings: fill_plot_options/data_options
      data_frame: get_data/data_frame
    out: [scatter_plot.json]
  show_scatter_plot:
    run: ui-outputs/plotly.cwl
    in:
      a: plot_scatter/scatter_plot.json
    out: []
    label: 'Display Scatter Plot'
