import React from 'react';
import {
  defaultContext,
  VulcanProvider
} from '../../../../contexts/vulcan_context';
import renderer from 'react-test-renderer';
import InputFeed from '../input_feed';
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
  setSession
} from '../../../../actions/vulcan_actions';
import Card from '@material-ui/core/Card';

describe('InputFeed', () => {
  const createNodeMock = (element: any) => {
    if (element.type === 'textarea') {
      return document.createElement('textarea');
    } else {
      return null;
    }
  };

  it('renders complete UI steps and error steps', () => {
    const workflow = createWorkflowFixture({
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
          })
        )
      ),
      setDownloadedData('https://download1', 'default-value'),
      setSession({
        project_name: 'test-project',
        workflow_name: workflow.name,
        key: 'test',
        inputs: {
          'second/response': 'default-value',
          'third/response': 'default-value'
        },
        reference_figure_id: 1
      })
    ]);

    const component = renderer.create(
      <VulcanProvider
        state={state}
        useActionInvoker={defaultContext.useActionInvoker}
      >
        <InputFeed />
      </VulcanProvider>,
      {createNodeMock}
    );

    let instance = component.root;

    expect(instance.findAllByProps({className: 'step-error'}).length).toEqual(
      1
    );

    expect(
      instance
        .findAllByType(Card)
        .filter((t) => t.props.className.endsWith('step-user-input')).length
    ).toEqual(2);

    expect(component.toJSON()).toMatchSnapshot();
  });
});
