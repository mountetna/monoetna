import React, {useContext} from 'react';

import {VulcanContext} from '../../../contexts/vulcan';
import {OUTPUT_COMPONENT} from '../../../steps';

import StepName from './step_name';
import RawOutput from '../user_interactions/outputs/raw';
import LinkOutput from '../user_interactions/outputs/link';
import PlotOutput from '../user_interactions/outputs/plot';
import ConsignmentOutput from '../user_interactions/outputs/consignment';

import {statusOfStep, uiOutputOfStep, uiStepInputDataLink, uiStepInputDataRaw} from "../../../selectors/workflow";
import {WorkflowStep} from "../../../api/types";

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

  let data;
  if (['default', OUTPUT_COMPONENT.LINK].includes(stepType)) {
    data = uiStepInputDataLink(step, state.status);
  } else {
    data = uiStepInputDataRaw(step, state.status, state.data);
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
