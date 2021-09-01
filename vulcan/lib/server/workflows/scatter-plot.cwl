cwlVersion: v1.1
class: Workflow

inputs:
  1_Filler__spacer:
    type: int
    default: 200
    label: "Filler - will be replaced with a data load step!"
    doc: "Does nothing"

outputs:
  the_plot:
    type: File
    outputSource: plot_scatter/plot.json

steps:
  get_data:
    run: scripts/mockDF.cwl
    label: 'Fetch Data'
    in:
      a: 1_Filler__spacer
    out: [data_frame]
  fill_plot_options:
    run: ui-queries/scatter-plotly.cwl
    label: 'Set plot options'
    in:
      data_frame: get_data/data_frame
    out: [plot_setup]
  plot_scatter:
    run: scripts/make_scatter.cwl
    label: 'Create Scatter Plot'
    in:
      plot_setup: fill_plot_options/plot_setup
      data_frame: get_data/data_frame
    out: [plot.json]
  show_scatter_plot:
    run: ui-outputs/plotly.cwl
    in:
      a: plot_scatter/plot.json
    out: []
    label: 'Display Scatter Plot'
