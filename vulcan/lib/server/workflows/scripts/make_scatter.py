from archimedes.functions.dataflow import output_path, input_path, json, input_json
from archimedes.functions.plotting import pio, scatter_plotnine
from archimedes.functions.utils import pandas as pd

df = pd.read_json(input_path("data_frame"), dtype = False)
plot_setup = input_json("plot_setup")

del_keyes = []
for key in plot_setup.keys():
    if plot_setup[key]=="make":
        del_keyes.extend([key])

if len(del_keyes)>0:
    for key in del_keyes:
        del plot_setup[key]

# Make & output plot
fig = scatter_plotnine(
    df, **plot_setup)

fig.save(output_path('plot.pdf'))