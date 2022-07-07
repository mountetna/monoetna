library(dataflow)

str <- input_str('string')
num <- input_num('number')
int <- input_int('integer')
bool <- input_bool('boolean')

output_var(str, 'str')
output_var(num, 'num')
output_var(int, 'int')
output_var(bool, 'bool')