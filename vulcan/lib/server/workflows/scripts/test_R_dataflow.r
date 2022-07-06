library(dataflow)

str <- input_str('string')
num <- input_num('number')
int <- input_str('integer')
bool <- input_str('boolean')

output_var(str, 'str')
output_var(num, 'num')
output_var(int, 'int')
output_var(bool, 'bool')