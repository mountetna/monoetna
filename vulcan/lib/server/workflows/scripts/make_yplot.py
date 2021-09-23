from archimedes.functions.dataflow import output_path, input_path, json, input_json
from archimedes.functions.plotting import pio, y_plotly
from archimedes.functions.utils import pandas as pd

df = pd.read_json(input_path("data_frame"), dtype = False)
plot_setup = input_json("plot_setup")

# Make & output plot
fig = y_plotly(
    df, **plot_setup)

with open(output_path('plot.json'), 'w') as output_file:
    json.dump(json.loads(pio.to_json(fig)), output_file)
