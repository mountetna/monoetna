import {
  FETCH_WORKFLOWS,
  FETCH_WORKFLOW,
  SUBMIT_INPUTS
} from '../../actions/archimedes_actions';
import reducer from '../archimedes_reducer';

describe('Archimedes reducer', () => {
  it('should handle SUBMIT_INPUTS and include data URLs and status', () => {
    expect(
      reducer(
        {
          workflow: {
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
        },
        {
          type: 'SUBMIT_INPUTS',
          status: {
            first_step: {
              status: 'complete',
              data_url: `${CONFIG.archimedes_host}/data1`
            },
            ui_pick_subset: {
              status: 'pending'
            },
            final_step: {
              status: 'pending',
              data_url: null
            }
          }
        }
      )
    ).toEqual({
      workflow: {
        first_step: {
          in: [],
          out: ['choice_set'],
          run: 'first_step.cwl',
          status: 'complete',
          data_url: `${CONFIG.archimedes_host}/data1`
        },
        ui_pick_subset: {
          in: {
            all_pool_names: 'first_step/choice_set'
          },
          out: ['subset'],
          run: 'ui_pick_subset.cwl',
          status: 'pending'
        },
        final_step: {
          in: {
            bool_input: 'bool_input',
            data: 'ui_pick_subset/subset',
            int_input: 'int_input'
          },
          out: ['sample_data'],
          run: 'final_step.cwl',
          status: 'pending',
          data_url: null
        }
      }
    });
  });
});
