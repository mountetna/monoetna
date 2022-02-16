from archimedes.functions.dataflow import output_path, input_path, json, input_json
from archimedes.functions.plotting import pio, scatter_plotly
from archimedes.functions.utils import pandas as pd
from archimedes.functions.ui_complements import subsetDF_index_targets

df = pd.read_json(input_path("data_frame"), dtype = False)
plot_setup = input_json("plot_setup")

if plot_setup['rows_use'] != {}:
    plot_setup['rows_use'] = subsetDF_index_targets(df, plot_setup['rows_use'])
else:
    plot_setup.pop('rows_use')

# Make & output plot
fig = scatter_plotly(
    df, **plot_setup)

with open(output_path('plot.json'), 'w') as output_file:
    json.dump(json.loads(pio.to_json(fig)), output_file)
