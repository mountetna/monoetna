import React from 'react';
import {VulcanProvider} from '../../../../contexts/vulcan';
import renderer from 'react-test-renderer';
import InputFeed from '../input_feed';

describe('InputFeed', () => {
  let state;

  beforeEach(() => {
    state = {
      workflow: {
        inputs: {},
        steps: [
          [
            {
              name: 'zero',
              run: 'scripts/something-wrong.cwl'
            },
            {
              name: 'first',
              run: 'scripts/something.cwl',
              out: ['output']
            },
            {
              name: 'second',
              run: 'ui-queries/ask-the-user.cwl',
              in: [{id: 'a', source: ['first', 'output']}],
              out: ['response']
            },
            {
              name: 'third',
              run: 'ui-queries/ask-again.cwl',
              in: [{id: 'a', source: ['first', 'output']}],
              out: ['response']
            },
            {
              name: 'fourth',
              run: 'ui-outputs/show-the-user.cwl'
            }
          ]
        ]
      },
      pathIndex: 0,
      session: {
        key: 'session_key',
        inputs: {
          'second/response': null,
          'third/response': null
        }
      },
      status: [
        [
          {
            name: 'zero',
            status: 'error',
            message: 'oops!'
          },
          {
            name: 'first',
            status: 'complete',
            downloads: {
              output: 'https://foo'
            },
            data: {
              output: [1, 2, 3]
            }
          },
          {
            name: 'second',
            status: 'pending'
          },
          {
            name: 'third',
            status: 'pending'
          },
          {
            name: 'fourth',
            status: 'pending'
          }
        ]
      ]
    };
  });

  it('renders complete UI steps and error steps', () => {
    // Wrap with Provider here so store gets passed down to child components in Context
    const component = renderer.create(
      <VulcanProvider state={state}>
        <InputFeed />
      </VulcanProvider>
    );

    let instance = component.root;

    expect(instance.findAllByProps({className: 'step-error'}).length).toEqual(
      1
    );

    expect(
      instance.findAllByProps({className: 'step-user-input step'}).length
    ).toEqual(2);

    expect(component.toJSON()).toMatchSnapshot();
  });
});
