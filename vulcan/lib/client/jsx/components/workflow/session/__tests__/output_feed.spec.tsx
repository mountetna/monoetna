import React from 'react';
import {Provider} from 'react-redux';
import {mockStore} from 'etna-js/spec/helpers';
import {VulcanProvider} from '../../../../contexts/vulcan_context';
import renderer from 'react-test-renderer';
import OutputFeed from '../output_feed';
import {
  createStatusFixture,
  createStepFixture,
  createStepStatusFixture,
  createWorkflowFixture
} from "../../../../test_utils/fixtures";
import {stateFromActions} from "../../../../test_utils/state";
import {setDownloadedData, setStatus, setWorkflow, setWorkflows} from "../../../../actions/vulcan";

describe('OutputFeed', () => {
  it('renders UI output steps', () => {
    // Link uses Redux store, so also include this Provider.
    const workflow = createWorkflowFixture({
      inputs: {}, steps: [
        [
          createStepFixture({
            name: 'first',
            run: 'scripts/something.cwl',
            out: ['output']
          }),
          createStepFixture({
            name: 'second',
            run: 'ui-outputs/plotly.cwl',
            in: [{id: 'a', source: 'first/output'}],
          }),
          createStepFixture({
            name: 'third',
            run: 'scripts/make-link.cwl',
            out: ['output']
          }),
          createStepFixture({
            name: 'fourth',
            run: 'ui-outputs/link.cwl',
            in: [{id: 'a', source: 'third/output'}],
          }),
          createStepFixture({
            name: 'fifth',
            run: 'scripts/make-raw.cwl',
            out: ['output']
          }),
          createStepFixture({
            name: 'sixth',
            run: 'ui-outputs/raw.cwl',
            in: [{id: 'a', source: 'fifth/output'}],
          }),
        ]
      ]
    });

    const {state} = stateFromActions([
      setWorkflows([workflow]),
      setWorkflow(workflow),
      setStatus(createStatusFixture(workflow,
          createStepStatusFixture({
            name: 'first',
            status: 'complete',
            downloads: {
              output: 'https://foo2'
            },
          }),
          createStepStatusFixture({
            name: 'second',
            status: 'complete'
          }),
          createStepStatusFixture({
            name: 'third',
            status: 'complete',
            downloads: {
              output: 'https://foo'
            }
          }),
          createStepStatusFixture({
            name: 'fourth',
            status: 'complete'
          }),
          createStepStatusFixture({
            name: 'fifth',
            status: 'complete',
            downloads: {
              output: 'https://foo'
            },
          }),
          createStepStatusFixture({
            name: 'sixth',
            status: 'complete'
          }),
      )),
      setDownloadedData('https://foo', ['Foo data']),
      setDownloadedData('https://foo2', ['Foo2 data'])
    ]);

    const store = mockStore(state);
    // Wrap with Provider here so store gets passed down to child components in Context
    const component = renderer.create(
        <Provider store={store}>
          <VulcanProvider state={state}>
            <OutputFeed />
        </VulcanProvider>
      </Provider>
    );

    let instance = component.root;

    expect(instance.findAllByProps({className: 'step-output'}).length).toEqual(
      3
    );

    expect(component.toJSON()).toMatchSnapshot();
  });
});
