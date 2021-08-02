import React, {useContext, useEffect, useMemo, useState} from 'react';
import { ComponentStory, ComponentMeta } from '@storybook/react';
import PrimaryInputs from "../lib/client/jsx/components/workflow/session/primary_inputs";
import {useWorkflowUtils, workflowUtilsBuilder} from "../lib/client/jsx/test_utils/workflow_utils";
import {defaultWorkflow, TYPE, Workflow, WorkflowInput, WorkflowStep} from "../lib/client/jsx/api_types";
import { DataEnvelope } from '../lib/client/jsx/components/workflow/user_interactions/inputs/input_types';
import StepUserInput from '../lib/client/jsx/components/workflow/steps/step_user_input';
import { useWorkflow } from '../lib/client/jsx/contexts/workflow_context';
import { VulcanContext } from '../lib/client/jsx/contexts/vulcan_context';
import {mapSome, Maybe, some, withDefault} from "../lib/client/jsx/selectors/maybe";
import {isPendingUiQuery} from "../lib/client/jsx/selectors/workflow_selectors";

interface Parameterization {
  cwlParams?: DataEnvelope<any>,
  type: string
}

function ParameterizedStepUserInput({cwlParams = {}, type}: Parameterization) {
  const utils = useWorkflowUtils();
  const [step, setState] = useState(null as Maybe<WorkflowStep>);
  const {state} = useContext(VulcanContext);
  const {status, data, session} = state;

  useEffect(() => {
    utils.setWorkflow('test');
    Object.keys(cwlParams).forEach(paramName => {
      utils.addStep(paramName, { out: ['output'] });
    });
    setState(some(utils.addStep('test-inputs', {
        run: `ui-queries/${type}`,
        in: Object.keys(cwlParams).map(paramName => ({source: `${paramName}/output`, id: paramName})),
        out: ['result'],
    })));

    Object.keys(cwlParams).forEach(paramName => {
      utils.forceDownloadedData(`${paramName}/output`, cwlParams[paramName]);
    });
  }, [utils, cwlParams, type]);

  return withDefault(mapSome(step, step => {
    if (!isPendingUiQuery(step, status, data, session)) {
      return <div>
        Not all data inputs are loaded.
      </div>;
    }

    return <StepUserInput step={step} hideLabel={false}/>;
  }), null);
}

export default {
  title: 'Inputs/StepUserInput',
  component: ParameterizedStepUserInput,
} as ComponentMeta<typeof ParameterizedStepUserInput>;

const Template: ComponentStory<typeof ParameterizedStepUserInput> = (args: any) => <ParameterizedStepUserInput {...args}/>;

export const ScatterPlotlyFULL = Template.bind({});
ScatterPlotlyFULL.args = {
  type: TYPE.SCATTER_PLOTLY,
  cwlParams: {
      'data_frame': require('./mockDF.json')
  }
};

export const ScatterPlotlyUMAP = Template.bind({});
ScatterPlotlyUMAP.args = {
  type: TYPE.SCATTER_PLOTLY,
  cwlParams: {
      'data_options': require('./color_options.json'),
      'hide': ['x_by', 'y_by']
  }
};
