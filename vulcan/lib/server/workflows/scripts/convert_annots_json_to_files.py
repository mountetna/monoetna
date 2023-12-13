from archimedes.functions.utils import pandas as pd
from archimedes.functions.dataflow import input_path, output_path

df = pd.read_json(input_path("annots.json"), dtype = False)
if df.shape[1] < 2 or not any([i.startswith('annot') for i in df.columns]):
    raise Exception("No annotation columns added. Return to previous step!")
df.to_excel(output_path("annots.xlsx"), sheet_name='sheet1', index=False)
df.to_csv(output_path("annots.csv"))
