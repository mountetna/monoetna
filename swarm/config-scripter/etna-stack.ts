import * as yaml from 'yaml';
import {prepYamlInputs} from './input-yaml';

const {obj, values: [name, env]} = prepYamlInputs(['', ''])

console.log(yaml.stringify(obj, {}))
