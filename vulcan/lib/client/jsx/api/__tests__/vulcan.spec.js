const fs = require('fs');
const path = require('path');

import nock from 'nock';

import {getWorkflows, getWorkflow, submitInputs, getData} from '../vulcan';
import {mockStore, stubUrl, mockFetch, cleanStubs} from 'etna-js/spec/helpers';

describe('Vulcan API', () => {
  afterEach(() => {
    cleanStubs();
    nock.cleanAll();
  });
  afterAll(nock.restore);

  mockFetch();

  it('getWorkflows returns all workflows', (done) => {
    const expectedWorkflows = {
      analysis_1: {},
      analysis_2: {},
      analysis_3: {}
    };

    stubUrl({
      verb: 'get',
      path: ROUTES.fetch_workflows(),
      response: expectedWorkflows,
      headers: {
        'Content-type': 'application/json'
      },
      host: CONFIG.vulcan_host
    });

    return getWorkflows().then((data) => {
      expect(data).toEqual(expectedWorkflows);
      done();
    });
  });

  it('getWorkflow returns a JSON workflow', (done) => {
    const sample_yaml = fs.readFileSync(
      path.resolve(__dirname, '../../spec/fixtures/sample_cwl.yaml'),
      'utf8'
    );

    // This should match the CWL file in `lib/client/jsx/spec/fixtures/sample_cwl.yaml`
    const expectedWorkflow = {
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
      steps: [
        {
          in: [],
          out: ['choice_set'],
          run: 'first_step.cwl',
          name: 'first_step'
        },
        {
          in: {
            all_pool_names: 'first_step/choice_set'
          },
          out: ['subset'],
          run: 'ui_pick_subset.cwl',
          name: 'ui_pick_subset'
        },
        {
          in: {
            bool_input: 'bool_input',
            data: 'ui_pick_subset/subset',
            int_input: 'int_input'
          },
          out: ['sample_data'],
          run: 'final_step.cwl',
          name: 'final_step'
        }
      ]
    };

    stubUrl({
      verb: 'get',
      path: ROUTES.fetch_workflow('test'),
      response: sample_yaml,
      headers: {
        'Content-type': 'text/yaml'
      },
      host: CONFIG.vulcan_host
    });

    return getWorkflow('test').then((data) => {
      expect(data).toEqual(expectedWorkflow);
      done();
    });
  });

  it('submitInputs correctly converts the input Array', (done) => {
    const inputs = [
      {
        name: 'first_step'
      },
      {
        all_pool_names: [1, 2, 3, 4],
        name: 'ui_pick_subset'
      },
      {
        name: 'final_step'
      }
    ];

    const status = {
      first_step: {
        status: 'complete',
        data_url: 'https://something'
      },
      ui_pick_subset: {
        status: 'pending'
      },
      final_step: {
        status: 'pending'
      }
    };

    stubUrl({
      verb: 'post',
      path: ROUTES.submit_inputs('test'),
      request: {
        first_step: {
          name: 'first_step'
        },
        ui_pick_subset: {
          all_pool_names: [1, 2, 3, 4],
          name: 'ui_pick_subset'
        },
        final_step: {
          name: 'final_step'
        },
        project_name: CONFIG.project_name
      },
      response: status,
      headers: {
        'Content-type': 'application/json'
      },
      host: CONFIG.vulcan_host
    });

    return submitInputs('test', inputs).then((data) => {
      expect(data).toEqual(status);
      done();
    });
  });

  it('getData returns data from the specified location', (done) => {
    const url = '/api/workflows/test/file1.txt';
    const data = [1, 2, 4, 'abc'];

    stubUrl({
      verb: 'get',
      path: url,
      response: data,
      headers: {
        'Content-type': 'application/json'
      },
      host: CONFIG.vulcan_host
    });

    return getData(url).then((returnedData) => {
      expect(data).toEqual(returnedData);
      done();
    });
  });
});
