from archimedes.functions.dataflow import output_path, input_path, json
from archimedes.functions.matrix import pd
from archimedes.functions.utils import random
from archimedes.functions.plotting import px, pio

a = open(input_path('a'), 'r').read()

NUM_POINTS = 1000000
NUM_FEATURES = 3


df = pd.DataFrame(random.randint(0, NUM_POINTS / 5,
                                 size=(NUM_POINTS, NUM_FEATURES)), columns=list('ABC'))

fig = px.scatter(df, x="A", y="B", color="C")

with open(output_path('umap'), 'w') as output_file:
    json.dump(json.loads(pio.to_json(fig)), output_file)
