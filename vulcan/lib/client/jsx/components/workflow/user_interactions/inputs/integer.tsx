import {WithInputParams} from './input_types';
import {NumberInput} from './float';

export default function IntegerInput({onChange, label, minWidth, data, ...props}: WithInputParams<{label?: string, minWidth?: number}, number | null, number | null>) {
  return NumberInput({onChange, label, minWidth, data, ...props}, true)
}
