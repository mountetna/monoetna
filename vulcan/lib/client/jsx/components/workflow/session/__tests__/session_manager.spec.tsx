import React from 'react';
import {VulcanProvider} from '../../../../contexts/vulcan_context';
import {rest} from 'msw';
import {setupServer} from 'msw/node';
import {render, fireEvent, waitFor, screen} from '@testing-library/react';
import '@testing-library/jest-dom/extend-expect';
import SessionManager from '../session_manager';
import {stateFromActions} from '../../../../test_utils/state';
import {
  createStatusFixture,
  createStepFixture,
  createStepStatusFixture,
  createWorkflowFixture
} from '../../../../test_utils/fixtures';
import {
  setDownloadedData,
  setStatus,
  setWorkflow,
  setWorkflows,
  addValidationErrors
} from '../../../../actions/vulcan';

describe('SessionManager', () => {
  const workflow = createWorkflowFixture({
    projects: ['test-project'],
    inputs: {},
    name: 'test.cwl',
    steps: [
      [
        createStepFixture({name: 'zero'}),
        createStepFixture({name: 'first', out: ['output']}),
        createStepFixture({name: 'second', out: ['output']}),
        createStepFixture({
          name: 'third',
          run: 'ui-queries/nested-select-autocomplete.cwl',
          in: [{id: 'a', source: 'first/output'}],
          out: ['response']
        }),
        createStepFixture({
          name: 'fourth',
          run: 'ui-queries/multiple-multiselect-string-all.cwl',
          in: [{id: 'a', source: 'second/output'}],
          out: ['response']
        }),
        createStepFixture({
          name: 'fifth',
          run: 'ui-outputs/show-the-user.cwl'
        })
      ]
    ]
  });

  const handlers = [
    rest.post(
      'https://vulcan.test/api/test-project/session/test.cwl/status',
      (req, res, ctx) => {
        return res(ctx.json({success: true}));
      }
    ),
    rest.get('https://vulcan.test/api/workflows', (req, res, ctx) => {
      return res(ctx.json([workflow]));
    })
  ];
  const server = setupServer(...handlers);

  const invoke = jest.fn();

  beforeAll(() => server.listen());
  afterEach(() => server.resetHandlers());
  afterAll(() => server.close());

  it('submits the session if no validation errors', async () => {
    server.use(
      rest.post(
        'https://vulcan.test/api/test-project/session/test.cwl',
        (req, res, ctx) => {
          return res(ctx.json({success: true}));
        }
      )
    );

    const {state} = stateFromActions([
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
      ),
      setDownloadedData('https://download1', {abc: null, '123': {xyz: null}}),
      setDownloadedData('https://download2', {
        abc: ['1', '2'],
        '123': ['a', 'b']
      })
    ]);

    render(
      <VulcanProvider state={state} useActionInvoker={() => invoke}>
        <SessionManager />
      </VulcanProvider>
    );

    await waitFor(() => screen.getByText('Run'));

    expect(screen.getByText('Run')).not.toBeDisabled();

    fireEvent.click(screen.getByText('Run'));

    expect(invoke).not.toHaveBeenCalled();
  });

  it('shows Validation errors on click Run button', async () => {
    const {state} = stateFromActions([
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
      ),
      setDownloadedData('https://download1', {abc: null, '123': {xyz: null}}),
      setDownloadedData('https://download2', {
        abc: ['1', '2'],
        '123': ['a', 'b']
      }),
      addValidationErrors('third', 'third', ['Oops!'])
    ]);

    render(
      <VulcanProvider state={state} useActionInvoker={() => invoke}>
        <SessionManager />
      </VulcanProvider>
    );

    await waitFor(() => screen.getByText('Run'));

    expect(screen.getByText('Run')).not.toBeDisabled();

    fireEvent.click(screen.getByText('Run'));

    expect(invoke).toHaveBeenCalled();
  });
});
