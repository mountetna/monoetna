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
    run: ui-outputs/scatter-plotly.cwl
    label: 'Set plot options'
    in:
      df: get_data/data_frame
      options: get_data/data_options
    out: [plot_settings]
  plot_scatter:
    run: ui-queries/scatter_make_scatter.cwl
    label: 'Create Scatter Plot'
    in:
      a: fill_plot_options/plot_settings
      df: get_data/data_frame
    out: [scatter_plot.json]
  show_scatter_plot:
    run: ui-outputs/plotly.cwl
    in:
      a: plot_scatter/scatter_plot.json
    out: []
    label: 'Display Scatter Plot'
