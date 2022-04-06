from archimedes.functions.dataflow import input_path
from archimedes.functions.utils import pandas as pd
from archimedes.functions.plotting import output_column_types

df = pd.read_json(input_path("data_frame"), dtype = False)

output_column_types(df, "continuous_cols", "discrete_cols")
