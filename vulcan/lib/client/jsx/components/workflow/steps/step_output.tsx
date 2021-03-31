import React, {useContext} from 'react';

import {VulcanContext} from '../../../contexts/vulcan_context';
import {OUTPUT_COMPONENT} from '../../../api_types';

import StepName from './step_name';
import RawOutput from '../user_interactions/outputs/raw';
import LinkOutput from '../user_interactions/outputs/link';
import PlotOutput from '../user_interactions/outputs/plot';
import ConsignmentOutput from '../user_interactions/outputs/consignment';

import {statusOfStep, uiOutputOfStep, stepInputDataUrls, stepInputDataRaw} from "../../../selectors/workflow_selectors";
import {WorkflowStep} from "../../../api_types";

const OUTPUTS = {
  default: LinkOutput,
  [OUTPUT_COMPONENT.LINK]: LinkOutput,
  [OUTPUT_COMPONENT.PLOTLY]: PlotOutput,
  [OUTPUT_COMPONENT.CONSIGNMENT]: ConsignmentOutput,
  [OUTPUT_COMPONENT.RAW]: RawOutput
};

export default function StepOutput({step}: {step: WorkflowStep}) {
  let {state} = useContext(VulcanContext);
  const stepStatus = statusOfStep(step, state.status);
  const uiOutput = uiOutputOfStep(step);

  if (!stepStatus || !uiOutput) return null;
  const stepType = uiOutput in OUTPUTS ? uiOutput : 'default';
  // console.log({stepType});

  let data;
  if (['default', OUTPUT_COMPONENT.LINK].includes(stepType)) {
    data = stepInputDataUrls(step, state.status);
  } else {
    data = stepInputDataRaw(step, state.status, state.data, state.session);
  }

  let OutputComponent = OUTPUTS[stepType];

  return (
    <div className='step-output'>
      <StepName step={step} />
      <div className='outputs-pane'>
        <OutputComponent data={data}/>
      </div>
    </div>
  );
}
