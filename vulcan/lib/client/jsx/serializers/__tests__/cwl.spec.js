const fs = require('fs');
const path = require('path');

import CwlSerializer from '../cwl';

describe('CwlSerializer', () => {
  it('converts a yaml doc to json', () => {
    const sample_yaml = fs.readFileSync(
      path.resolve(__dirname, '../../spec/fixtures/sample_cwl.yaml'),
      'utf8'
    );
    const cwl = new CwlSerializer(sample_yaml);

    const doc = cwl.json;

    expect(Object.keys(doc)).toEqual([
      'cwlVersion',
      'class',
      'inputs',
      'outputs',
      'steps'
    ]);
    expect(Object.keys(doc.steps)).toEqual([
      'first_step',
      'ui_pick_subset',
      'final_step'
    ]);
  });
});
