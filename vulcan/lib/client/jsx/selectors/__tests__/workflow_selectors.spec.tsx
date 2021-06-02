import {isPendingUiQuery} from '../workflow_selectors';
import {VulcanState} from '../../reducers/vulcan_reducer';
import {WorkflowStep, Workflow} from '../../api_types';
import {
  createStatusFixture,
  createStepFixture,
  createStepStatusFixture,
  createWorkflowFixture
} from '../../test_utils/fixtures';
import {
  setDownloadedData,
  setStatus,
  setWorkflow,
  setWorkflows
} from '../../actions/vulcan';
import {stateFromActions} from '../../test_utils/state';

describe('Workflow Selectors', () => {
  let status: VulcanState['status'];
  let data: VulcanState['data'];
  let session: VulcanState['session'];
  let step: WorkflowStep;
  let workflow: Workflow;
  let state: VulcanState;

  describe('isPendingUiQuery', () => {
    beforeEach(() => {
      workflow = createWorkflowFixture({
        projects: ['test-project'],
        inputs: {},
        steps: [
          [
            createStepFixture({name: 'zero'}),
            createStepFixture({name: 'first', out: ['output']}),
            createStepFixture({
              name: 'second',
              run: 'ui-queries/ask-the-user.cwl',
              in: [{id: 'a', source: 'first/output'}],
              out: ['response']
            }),
            createStepFixture({
              name: 'third',
              run: 'ui-queries/ask-again.cwl',
              in: [{id: 'a', source: 'first/output'}],
              out: ['response']
            }),
            createStepFixture({
              name: 'fourth',
              run: 'ui-outputs/show-the-user.cwl'
            })
          ]
        ]
      });

      state = stateFromActions([
        setWorkflows([workflow]),
        setWorkflow(workflow, 'test-project'),
        setStatus(
          createStatusFixture(
            workflow,
            createStepStatusFixture({
              name: 'zero',
              status: 'error',
              error: 'Ooops!'
            }),
            createStepStatusFixture({
              name: 'first',
              status: 'complete',
              downloads: {output: 'https://download1'}
            })
          )
        ),
        setDownloadedData('https://download1', {abc: null, '123': null})
      ]);
    });

    it('returns true for ui-query with simple list input', () => {});

    it('return true for ui-query with nested hash input', () => {});

    it('returns false if data dependency not available yet', () => {});

    it('returns false if status is not pending', () => {});
  });
});
