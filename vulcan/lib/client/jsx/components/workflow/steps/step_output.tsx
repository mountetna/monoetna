import React, {useContext} from 'react';

import {VulcanContext} from '../../../contexts/vulcan_context';
import {OUTPUT_COMPONENT} from '../../../api_types';

import StepName from './step_name';
import RawOutput from '../user_interactions/outputs/raw';
import LinkOutput from '../user_interactions/outputs/link';
import { PlotlyOutput, PlotOutput, PngOutput } from '../user_interactions/outputs/plot';
import ConsignmentOutput from '../user_interactions/outputs/consignment';

import {statusOfStep, uiOutputOfStep, stepInputDataUrls, stepInputDataRaw} from "../../../selectors/workflow_selectors";
import {WorkflowStep} from "../../../api_types";

const OUTPUTS = {
  default: LinkOutput,
  [OUTPUT_COMPONENT.LINK]: LinkOutput,
  [OUTPUT_COMPONENT.PLOTLY]: PlotlyOutput,
  [OUTPUT_COMPONENT.PLOT]: PlotOutput,
  [OUTPUT_COMPONENT.PNG]: PngOutput,
  [OUTPUT_COMPONENT.CONSIGNMENT]: ConsignmentOutput,
  [OUTPUT_COMPONENT.RAW]: RawOutput
};

export default function StepOutput({step}: {step: WorkflowStep}) {
  let {state} = useContext(VulcanContext);
  const stepStatus = statusOfStep(step, state.status);
  const uiOutput = uiOutputOfStep(step);

  if (!stepStatus || !uiOutput) return null;
  const stepType = uiOutput in OUTPUTS ? uiOutput : 'default';

  let data;
  if (['default', OUTPUT_COMPONENT.LINK,OUTPUT_COMPONENT.PNG].includes(stepType)) {
    data = stepInputDataUrls(step, state.status);
  } else if ([OUTPUT_COMPONENT.PLOT].includes(stepType)){
    data = {
      url: stepInputDataUrls(step, state.status),
      raw: stepInputDataRaw(step, state.status, state.data, state.session)};
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
