from archimedes.functions.dataflow import output_path, input_tsv, json, input_json
from archimedes.functions.plotting import pio, scatter_plotly
from archimedes.functions.utils import pandas as pd

# Read inputs
df = input_tsv('data_frame')
plot_settings = input_json('plot_settings')

### Parse plot_settings inputs
# Change any "make" to None
def str_val_to_None(input: any, toNone: str = "make"):
    if (isinstance(input, str) and (input == toNone)):
        return None
    return input

plot_settings = dict(
    [
        key,
        str_val_to_None(plot_settings[key], "make")
    ] for key in plot_settings.keys())

# Make & output plot
fig = scatter_plotly(
    df, **plot_settings)

with open(output_path('scatter_plot.json'), 'w') as output_file:
    json.dump(json.loads(pio.to_json(fig)), output_file)