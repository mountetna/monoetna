import reducer from '../vulcan';
import {
  SET_DATA,
  SET_WORKFLOW,
  SET_WORKFLOWS,
  SET_STATUS
} from '../../actions/vulcan';

describe('Vulcan Reducer', () => {
  it('correctly sets all workflows', () => {
    const expectedWorkflows = {
      analysis_1: {},
      analysis_2: {},
      analysis_3: {}
    };

    const state = reducer(
      {},
      {
        type: SET_WORKFLOWS,
        workflows: expectedWorkflows
      }
    );

    expect(state).toEqual({
      workflows: expectedWorkflows
    });
  });

  it('correctly sets current workflow', () => {
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

    const state = reducer(
      {},
      {
        type: SET_WORKFLOW,
        workflow: expectedWorkflow
      }
    );

    expect(state).toEqual({workflow: expectedWorkflow});
  });

  it('correctly injects status updates', () => {
    const status = [
      [
        {
          name: 'first_step',
          status: 'complete',
          data_url: `${CONFIG.vulcan_host}/data/blobs/here`
        },
        {
          status: 'pending',
          data_url: null,
          name: 'ui_pick_subset'
        },
        {
          name: 'final_step',
          status: 'pending',
          data_url: null
        }
      ]
    ];

    const workflow = {
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
        [
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
      ]
    };

    const state = reducer({workflow}, {type: SET_STATUS, status});

    expect(state.workflow.steps).toEqual([
      [
        {
          in: [],
          out: ['choice_set'],
          run: 'first_step.cwl',
          name: 'first_step',
          status: 'complete',
          data_url: `${CONFIG.vulcan_host}/data/blobs/here`
        },
        {
          in: {
            all_pool_names: 'first_step/choice_set'
          },
          out: ['subset'],
          run: 'ui_pick_subset.cwl',
          name: 'ui_pick_subset',
          status: 'pending',
          data_url: null
        },
        {
          in: {
            bool_input: 'bool_input',
            data: 'ui_pick_subset/subset',
            int_input: 'int_input'
          },
          out: ['sample_data'],
          run: 'final_step.cwl',
          name: 'final_step',
          status: 'pending',
          data_url: null
        }
      ]
    ]);
  });

  it('correctly injects data payload to right step', () => {
    const workflow = {
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
        [
          {
            in: [],
            out: ['choice_set'],
            run: 'first_step.cwl',
            name: 'first_step',
            data_url: '/api/workflows/test/file1.txt'
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
      ]
    };

    const url = '/api/workflows/test/file1.txt';
    const data = [1, 2, 4, 'abc'];

    const state = reducer({workflow}, {type: SET_DATA, url, data});

    expect(state.workflow.steps).toEqual([
      [
        {
          in: [],
          out: ['choice_set'],
          run: 'first_step.cwl',
          name: 'first_step',
          data_url: '/api/workflows/test/file1.txt',
          data: [1, 2, 4, 'abc']
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
    ]);
  });

  it('correctly injects inputs into the session', () => {
    const session = {
      key: '123',
      project_name: CONFIG.project_name,
      workflow_name: 'test',
      inputs: []
    };

    const inputs = [[{a: 123, b: 321}]];

    const state = reducer({session}, {type: SET_INPUTS, inputs});

    expect(state.session.inputs).toEqual(inputs);
  });
});
