import React, {useContext} from 'react';

import {VulcanContext} from '../../../contexts/vulcan';
import {OUTPUT_COMPONENT} from '../../../models/steps';

import StepName from './step_name';
import RawOutput from '../user_interactions/outputs/raw';
import LinkOutput from '../user_interactions/outputs/link';
import PlotOutput from '../user_interactions/outputs/plot';
import ConsignmentOutput from '../user_interactions/outputs/consignment';

import {
  validStep,
  hasUiOutput,
  uiStepInputDataRaw,
  uiStepInputDataLink
} from '../../../utils/workflow';

export default function StepOutput({step, stepIndex}) {
  const {workflow, pathIndex, session, status} = useContext(VulcanContext);

  if (
    !validStep({workflow, pathIndex, stepIndex}) ||
    !session ||
    !status ||
    !step ||
    null === stepIndex ||
    !hasUiOutput(step)
  )
    return null;

  let stepType = step.run.split('/')[1].replace('.cwl', '');

  stepType = Object.values(OUTPUT_COMPONENT).includes(stepType)
    ? stepType
    : 'default';

  const OUTPUTS = {
    default: LinkOutput,
    [OUTPUT_COMPONENT.LINK]: LinkOutput,
    [OUTPUT_COMPONENT.PLOTLY]: PlotOutput,
    [OUTPUT_COMPONENT.CONSIGNMENT]: ConsignmentOutput,
    [OUTPUT_COMPONENT.RAW]: RawOutput
  };

  let data;

  if (['default', OUTPUT_COMPONENT.LINK].includes(stepType)) {
    data = uiStepInputDataLink({step, pathIndex, status});
  } else {
    data = uiStepInputDataRaw({step, pathIndex, status});
  }

  let OutputComponent = OUTPUTS[stepType];

  return (
    <div className='step-output'>
      <StepName
        step={step}
        status={status[pathIndex][stepIndex].status}
      ></StepName>
      <div className='outputs-pane'>
        <OutputComponent data={data}></OutputComponent>
      </div>
    </div>
  );
}
