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
    outputSource: make_plot/plot.json

steps:
  get_data:
    run: scripts/mockDF.cwl
    label: 'Fetch Data'
    in:
      a: 1_ignore__spacer
    out: [data_frame]
  fill_plot_options:
    run: ui-queries/bar-plotly.cwl
    label: 'Set plot options'
    in:
      data_frame: get_data/data_frame
    out: [plot_setup]
  make_plot:
    run: scripts/make_barplot.cwl
    label: 'Create Bar Plot'
    in:
      plot_setup: fill_plot_options/plot_setup
      data_frame: get_data/data_frame
    out: [plot.json]
  show_plot:
    run: ui-outputs/plotly.cwl
    in:
      a: make_plot/plot.json
    out: []
    label: 'Display Scatter Plot'
