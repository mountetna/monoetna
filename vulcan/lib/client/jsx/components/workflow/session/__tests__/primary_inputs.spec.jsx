import React from 'react';
import {VulcanProvider} from '../../../../contexts/vulcan';
import renderer from 'react-test-renderer';
import PrimaryInputs from '../primary_inputs';

describe('PrimaryInputs', () => {
  let state;

  beforeEach(() => {
    state = {
      workflow: {
        inputs: {
          anInt: {
            type: 'int',
            default: 1
          },
          aFloat: {
            type: 'float',
            default: 1.2
          },
          aBool: {
            type: 'boolean',
            default: true
          },
          anIntWithoutDefault: {
            type: 'int',
            default: null
          },
          aFloatWithoutDefault: {
            type: 'float',
            default: null
          },
          aBoolWithoutDefault: {
            type: 'boolean',
            default: null
          }
        }
      },
      session: {
        key: 'session_key'
      }
    };
  });

  it('renders each type of input', () => {
    // Wrap with Provider here so store gets passed down to child components in Context
    const component = renderer.create(
      <VulcanProvider state={state}>
        <PrimaryInputs />
      </VulcanProvider>
    );

    let instance = component.root;

    expect(instance.findAllByType('input').length).toEqual(
      Object.keys(state.workflow.inputs).length
    );

    expect(component.toJSON()).toMatchSnapshot();
  });
});
