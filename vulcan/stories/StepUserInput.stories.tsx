import React, {useContext, useEffect, useMemo} from 'react';
import { ComponentStory, ComponentMeta } from '@storybook/react';
import PrimaryInputs from "../lib/client/jsx/components/workflow/session/primary_inputs";
import {useWorkflowUtils, workflowUtilsBuilder} from "../lib/client/jsx/test_utils/workflow_utils";
import {defaultWorkflow, TYPE, WorkflowInput, WorkflowStep} from "../lib/client/jsx/api_types";
import { DataEnvelope } from '../lib/client/jsx/components/workflow/user_interactions/inputs/input_types';
import StepUserInput from '../lib/client/jsx/components/workflow/steps/step_user_input';
import { useWorkflow } from '../lib/client/jsx/contexts/workflow_context';
import { VulcanContext } from '../lib/client/jsx/contexts/vulcan_context';

interface Parameterization {
  cwlParams?: DataEnvelope<any>,
  type: string
}

function ParameterizedStepUserInput({cwlParams = {}, type}: Parameterization) {
  const utils = useMemo(() => workflowUtilsBuilder(), []);
  const utils2 = useWorkflowUtils();
  const {workflow: stateWorkflow} = useWorkflow();
  const {state, stateRef} = useContext(VulcanContext);

  const step = useMemo(() => {
    utils.setWorkflow('test');
    Object.keys(cwlParams).forEach(paramName => {
        utils.addStep(paramName, { out: ['output'] });
    });
    return utils.addStep('test-inputs', {
        run: `ui-queries/${type}`,
        in: Object.keys(cwlParams).map(paramName => ({source: `${paramName}/output`, id: paramName})),
        out: ['result'],
    });
  }, [utils, cwlParams]);

  useEffect(() => {
    utils2.setWorkflow('test', utils.workflow);
    let status = stateRef.current.status;
    Object.keys(cwlParams).forEach(paramName => {
        status = utils2.forceDownloadedData(`${paramName}/output`, cwlParams[paramName], status);
    });
  }, [utils, utils.workflow, cwlParams])

  if (stateWorkflow.name !== utils.workflow.name) return null;

  return <StepUserInput step={step} hideLabel={false}/>;
}

export default {
  title: 'Inputs/StepUserInput',
  component: ParameterizedStepUserInput,
} as ComponentMeta<typeof ParameterizedStepUserInput>;

const Template: ComponentStory<typeof ParameterizedStepUserInput> = (args: any) => <ParameterizedStepUserInput {...args}/>;

export const ScatterPlotly = Template.bind({});
ScatterPlotly.args = {
  type: TYPE.SCATTER_PLOTLY,
  cwlParams: {
      'data_frame': require('./mockDF.json'),
      'data_options': require('./color_options.json'),
  }
};
