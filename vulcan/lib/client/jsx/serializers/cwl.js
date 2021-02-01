import YAML from 'yaml';

import Serializer from './interface';

export default class CwlSerializer extends Serializer {
  get json() {
    // Convert a CWL YAML file into JSON
    const doc = YAML.parseDocument(this.raw);
    return doc.toJSON();
  }
}
