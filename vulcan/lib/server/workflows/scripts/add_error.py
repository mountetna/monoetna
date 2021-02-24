from archimedes.functions.dataflow import output_path, input_path

a = int(open(input_path('a'), 'r').read())
b = int(open(input_path('b'), 'r').read())

raise Exception(f"Something bad happened to this step! Inputs: ${a} + ${b}")
