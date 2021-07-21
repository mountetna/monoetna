import React, {useContext} from 'react';
import {defaultContext, VulcanContext, VulcanProvider} from '../../../../contexts/vulcan_context';
import renderer, {act} from 'react-test-renderer';
import PrimaryInputs from '../primary_inputs';
import {createWorkflowFixture} from "../../../../test_utils/fixtures";
import {awaitBefore, integrateElement, setupBefore} from "../../../../test_utils/integration";
import {useWorkflowUtils} from "../../../../test_utils/workflow_utils";

describe('PrimaryInputs', () => {
  const integrated = setupBefore(() => integrateElement(<PrimaryInputs/>));
  const workflowHelpers = setupBefore(() => integrated.value.runHook(() => useWorkflowUtils()));
  const contextData = setupBefore(() => integrated.value.runHook(() => useContext(VulcanContext)));

  awaitBefore(async () => {
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
          default: false,
          doc: 'help',
        }
      }
    });

    workflowHelpers.value.setWorkflow(workflow.name, workflow)
  })

  it('renders each type of input', async () => {
    const {stateRef} = contextData.value;
    const {node} = integrated.value;
    expect(node.root.findAllByType('input').length).toEqual(
        Object.keys(stateRef.current.workflow?.inputs || {}).length
    );

    expect(node.toJSON()).toMatchSnapshot();
  });
});
