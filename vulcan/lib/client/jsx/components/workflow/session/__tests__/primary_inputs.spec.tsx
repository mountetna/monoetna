import React from 'react';
import {defaultContext, VulcanProvider} from '../../../../contexts/vulcan_context';
import renderer from 'react-test-renderer';
import PrimaryInputs from '../primary_inputs';
import {stateFromActions} from "../../../../test_utils/state";
import {createWorkflowFixture} from "../../../../test_utils/fixtures";
import {setWorkflow, setWorkflows} from "../../../../actions/vulcan";

describe('PrimaryInputs', () => {
  it('renders each type of input', () => {
    const workflow = createWorkflowFixture({
      inputs: {
        anInt: {
          type: 'int',
          default: 1,
          doc: 'help',
        },
        aFloat: {
          type: 'float',
          default: 1.2,
          doc: 'help',
        },
        aBool: {
          type: 'boolean',
          default: true,
          doc: 'help',
        },
        anIntWithoutDefault: {
          type: 'int',
          default: null,
          doc: 'help',
        },
        aFloatWithoutDefault: {
          type: 'float',
          default: null,
          doc: 'help',
        },
        aBoolWithoutDefault: {
          type: 'boolean',
          default: null,
          doc: 'help',
        }
      }
    });

    const { state } = stateFromActions([
        setWorkflows([workflow]),
        setWorkflow(workflow, 'test'),
    ])

    const component = renderer.create(
        <VulcanProvider state={state} useActionInvoker={defaultContext.useActionInvoker}>
          <PrimaryInputs/>
        </VulcanProvider>
    );

    let instance = component.root;

    expect(instance.findAllByType('input').length).toEqual(
        Object.keys(workflow.inputs).length
    );

    expect(component.toJSON()).toMatchSnapshot();
  });
});
