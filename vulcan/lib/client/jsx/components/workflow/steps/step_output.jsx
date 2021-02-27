import React, {useContext} from 'react';
import Plot from 'react-plotly.js';

import {VulcanContext} from '../../../contexts/vulcan';

import StepName from './step_name';

import {
  validStep,
  hasUiOutput,
  uiStepInputDataRaw
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

  // We need to extract the data from the input source.
  let rawInputData = uiStepInputDataRaw({step, pathIndex, status});

  let Component;
  switch (step.run.split('/')[1].replace('.cwl', '')) {
    case 'plotly':
      // Plotly.js data payload should be in format of:
      // {
      //   data: <JSON>,
      //   layout: <JSON>
      // }
      if (null === rawInputData || !rawInputData.data || !rawInputData.layout)
        return null;
      Component = (
        <Plot data={rawInputData.data} layout={rawInputData.layout}></Plot>
      );
      break;
    case 'raw':
    default:
      break;
  }

  if (!Component) return null;

  return (
    <div className='step-output-plotly'>
      <StepName
        step={step}
        status={status[pathIndex][stepIndex].status}
      ></StepName>
      <div className='outputs-pane'>{Component}</div>
    </div>
  );
}
