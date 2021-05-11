from archimedes.functions.dataflow import output_path, input_path

a = int(open(input_path('a'), 'r').read())
b = int(open(input_path('b'), 'r').read())

with open(output_path('sum'), 'w') as output_file:
    output_file.write(str(a + b))
