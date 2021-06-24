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
} from '../../actions/vulcan_actions';
import {stateFromActions} from '../../test_utils/state';

describe('Workflow Selectors', () => {
  let workflow: Workflow;

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
    });

    it('returns true for ui-query with simple list input', () => {
      let state = stateFromActions([
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
        setDownloadedData('https://download1', ['1', '2', '3'])
      ])['state'];

      expect(
        isPendingUiQuery(
          workflow.steps[0][2],
          state.status,
          state.data,
          state.session
        )
      ).toEqual(true);
      expect(
        isPendingUiQuery(
          workflow.steps[0][3],
          state.status,
          state.data,
          state.session
        )
      ).toEqual(true);
    });

    it('return true for ui-query with nested hash input', () => {
      let state = stateFromActions([
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
        setDownloadedData('https://download1', {abc: null, '123': {xyz: null}})
      ])['state'];

      expect(
        isPendingUiQuery(
          workflow.steps[0][2],
          state.status,
          state.data,
          state.session
        )
      ).toEqual(true);
      expect(
        isPendingUiQuery(
          workflow.steps[0][3],
          state.status,
          state.data,
          state.session
        )
      ).toEqual(true);
    });

    it('returns false if data dependency not available yet', () => {
      let state = stateFromActions([
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
        )
      ])['state'];

      expect(
        isPendingUiQuery(
          workflow.steps[0][2],
          state.status,
          state.data,
          state.session
        )
      ).toEqual(false);
      expect(
        isPendingUiQuery(
          workflow.steps[0][3],
          state.status,
          state.data,
          state.session
        )
      ).toEqual(false);
    });

    it('returns false if status is not pending', () => {
      let state = stateFromActions([
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
            }),
            createStepStatusFixture({
              name: 'second',
              status: 'complete',
              downloads: {output: 'https://download2'}
            })
          )
        )
      ])['state'];

      expect(
        isPendingUiQuery(
          workflow.steps[0][2],
          state.status,
          state.data,
          state.session
        )
      ).toEqual(false);
    });
  });
});
