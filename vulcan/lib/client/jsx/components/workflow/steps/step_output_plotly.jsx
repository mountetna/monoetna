import React, {useContext} from 'react';
import Plot from 'react-plotly.js';

import {VulcanContext} from '../../../contexts/vulcan';

import StepName from './step_name';

import {
  validStep,
  hasUiOutput,
  uiStepInputDataRaw
} from '../../../utils/workflow';

export default function StepOutputPlotly({step, stepIndex}) {
  const {workflow, pathIndex, session, status, setInputs} = useContext(
    VulcanContext
  );

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

  // Plotly.js data payload should be in format of:
  // {
  //   data: <JSON>,
  //   layout: <JSON>
  // }
  if (null === rawInputData || !rawInputData.data || !rawInputData.layout)
    return null;

  return (
    <div className='step-output-plotly'>
      <StepName
        step={step}
        status={status[pathIndex][stepIndex].status}
      ></StepName>
      <div className='outputs-pane'>
        <Plot data={rawInputData.data} layout={rawInputData.layout}></Plot>
      </div>
      <div className=''>
        <p>Here you'll be able to download raw data</p>
      </div>
    </div>
  );
}
