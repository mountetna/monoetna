from archimedes.functions.dataflow import output_path, input_tsv, json, input_json
from archimedes.functions.plotting import pio, scatter_plotly
from archimedes.functions.utils import pandas as pd

# Read inputs
df = input_tsv('data_frame')
plot_settings = input_json('plot_settings')

### Parse plot_settings inputs
# Change any "make" to None
def str_val_to_None(input: any, toNone: str = "make"):
    if isinstance(input, str) and input == toNone:
        return None
    return input
plot_settings = dict(
    [
        key,
        str_val_to_None(plot_settings[key], "make")
    ] for key in plot_settings.keys())
# Sort inputs by how they need to be supplied
def move_keys(keys: list, orig: dict, out: dict = {}):
    for key in keys:
        if key in orig.keys():
            out[key] = orig[key]
            orig.pop(key)
    return out
direct_args = {
    'x_by': '0',
    'y_by': '1',
    'color_by': 'leiden',
    'color_order': 'increasing',
    'order_when_continuous_color': True,
    'size': 5
}
direct_args = move_keys(
    ['plot_title', 'legend_title', 'xlab', 'ylab'],
    plot_settings,
    direct_args
)
print(plot_settings)
print(direct_args)

px_args = plot_settings

# Make & output plot
fig = scatter_plotly(
    df, **direct_args,
    px_args = px_args)

with open(output_path('scatter_plot.json'), 'w') as output_file:
    json.dump(json.loads(pio.to_json(fig)), output_file)