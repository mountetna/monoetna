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
    expect(doc.steps.map((s) => s.name)).toEqual([
      'first_step',
      'ui_pick_subset',
      'final_step'
    ]);
  });

  it('handles steps with multiple dependencies', () => {
    // We'll feed the orderSteps method specific JSON, so don't
    //   need a real YAML input for the constructor.
    const cwl = new CwlSerializer('');

    const inputJson = {
      steps: {
        last_step: {
          in: {
            data1: 'second_step/output',
            data2: 'third_step/output'
          }
        },
        second_step: {
          in: {
            data: 'first_step/output'
          }
        },
        third_step: {
          in: {
            data: 'first_step/output2'
          }
        },
        first_step: {
          in: []
        }
      }
    };

    const output = cwl.orderSteps(inputJson);

    expect(output).toEqual({
      steps: [
        {
          in: [],
          name: 'first_step'
        },
        {
          in: {
            data: 'first_step/output'
          },
          name: 'second_step'
        },
        {
          in: {
            data: 'first_step/output2'
          },
          name: 'third_step'
        },
        {
          in: {
            data1: 'second_step/output',
            data2: 'third_step/output'
          },
          name: 'last_step'
        }
      ]
    });
  });
});
