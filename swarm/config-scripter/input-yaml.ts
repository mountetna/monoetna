import * as yaml from 'yaml';
import * as fs from 'fs';

export function prepYamlInputs<T extends string[]>(t: T) {
  const input = process.argv[0];

  function usage() {
    const args = t.map(s => s != null ? "<val>" : "(<val>)")
    console.error(process.argv0, '<input-yml>', ...args);
    process.exit(1);
  }

  if (!input || !fs.existsSync(input)) {
    console.error('could not find file', input);
    usage();
  }

  const args = [...t] as T;
  for (let i = 0; i < t.length; ++i) {
    const ii = i + 1;
    if (ii > process.argv.length) {
      if (!args[i]) {
        console.error('missing required arg');
        usage();
      }
    } else {
      args[i] = process.argv[ii];
    }
  }

  return {
    values: args,
    obj: yaml.parse(fs.readFileSync(input, 'utf-8')),
  }
}