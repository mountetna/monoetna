import React from 'react';
import {Provider} from 'react-redux';
import {mockStore} from 'etna-js/spec/helpers';
import {VulcanProvider} from '../../../../contexts/vulcan';
import renderer from 'react-test-renderer';
import OutputFeed from '../output_feed';

describe('OutputFeed', () => {
  let state;

  beforeEach(() => {
    state = {
      workflow: {
        inputs: {},
        steps: [
          [
            {
              name: 'first',
              run: 'scripts/something.cwl',
              out: ['output']
            },
            {
              name: 'second',
              run: 'ui-outputs/plotly.cwl',
              in: [{id: 'a', source: ['first', 'output']}],
              out: []
            },
            {
              name: 'third',
              run: 'scripts/make-link.cwl',
              in: [],
              out: ['output']
            },
            {
              name: 'fourth',
              run: 'ui-outputs/link.cwl',
              in: [{id: 'a', source: ['third', 'output']}],
              out: []
            },
            {
              name: 'fifth',
              run: 'scripts/make-raw.cwl',
              in: [],
              out: ['output']
            },
            {
              name: 'sixth',
              run: 'ui-outputs/raw.cwl',
              in: [{id: 'a', source: ['fifth', 'output']}],
              out: []
            }
          ]
        ]
      },
      pathIndex: 0,
      session: {
        key: 'session_key',
        inputs: {}
      },
      status: [
        [
          {
            name: 'first',
            status: 'complete',
            downloads: {
              output: 'https://foo'
            },
            data: {
              output: {data: [], layout: {}}
            }
          },
          {
            name: 'second',
            status: 'complete'
          },
          {
            name: 'third',
            status: 'complete',
            downloads: {
              output: 'https://foo'
            }
          },
          {
            name: 'fourth',
            status: 'complete'
          },
          {
            name: 'fifth',
            status: 'complete',
            downloads: {
              output: 'https://foo'
            },
            data: {
              output: 'some text!'
            }
          },
          {
            name: 'sixth',
            status: 'complete'
          }
        ]
      ]
    };
  });

  it('renders UI output steps', () => {
    // Link uses Redux store, so also include this Provider.
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
