const fs = require('fs');
const path = require('path');

import nock from 'nock';

import {fetchWorkflow} from '../archimedes_actions';
import {mockStore, stubUrl, mockFetch, cleanStubs} from 'etna-js/spec/helpers';

describe('Archimedes Actions', () => {
  afterEach(() => {
    cleanStubs();
    nock.cleanAll();
  });
  afterAll(nock.restore);

  mockFetch();

  it('fetchWorkflow dispatches a JSON object to the store', (done) => {
    const store = mockStore({});

    const sample_yaml = fs.readFileSync(
      path.resolve(__dirname, '../../spec/fixtures/sample_cwl.yaml'),
      'utf8'
    );

    // This should match the CWL file in `lib/client/jsx/spec/fixtures/sample_cwl.yaml`
    const expectedAction = {
      type: 'FETCH_WORKFLOW',
      workflow: {
        class: 'Workflow',
        cwlVersion: 'v1.1',
        inputs: {
          bool_input: {
            default: true,
            label: 'Sample boolean',
            type: 'boolean'
          },
          int_input: {
            default: 42,
            label: 'Sample input',
            type: 'int'
          }
        },
        outputs: {
          sample_data: {
            outputSource: 'final_step/sample_data',
            type: 'File'
          }
        },
        steps: {
          first_step: {
            in: [],
            out: ['choice_set'],
            run: 'first_step.cwl'
          },
          ui_pick_subset: {
            in: {
              all_pool_names: 'first_step/choice_set'
            },
            out: ['subset'],
            run: 'ui_pick_subset.cwl'
          },
          final_step: {
            in: {
              bool_input: 'bool_input',
              data: 'ui_pick_subset/subset',
              int_input: 'int_input'
            },
            out: ['sample_data'],
            run: 'final_step.cwl'
          }
        }
      }
    };

    stubUrl({
      verb: 'get',
      path: '/api/workflows/test',
      response: sample_yaml,
      headers: {
        'Content-type': 'text/yaml'
      },
      host: CONFIG.archimedes_host
    });

    return store.dispatch(fetchWorkflow('test')).then(() => {
      const actions = store.getActions();

      expect(actions[actions.length - 1]).toEqual(expectedAction);
      done();
    });
  });
});
