from archimedes.functions.dataflow import output_path, input_path, json, input_json
from archimedes.functions.plotting import pio, scatter_plotly, y_plotly, bar_plotly
from archimedes.functions.utils import pandas as pd
from archimedes.functions.ui_complements import subsetDF_index_targets

df = pd.read_json(input_path("data_frame"), dtype = False)
plot_setup = input_json("plot_setup")

if plot_setup['rows_use'] != {}:
    plot_setup['rows_use'] = subsetDF_index_targets(df, plot_setup['rows_use'])
else:
    plot_setup.pop('rows_use')

viz_fxn = {
    'scatter_plot': scatter_plotly,
    'y_plot': y_plotly,
    'bar_plot': bar_plotly
}[plot_setup['plot_type']]
plot_setup.pop('plot_type')

# Make & output plot
fig = viz_fxn(
    df, **plot_setup)

with open(output_path('plot.json'), 'w') as output_file:
    json.dump(json.loads(pio.to_json(fig)), output_file)
