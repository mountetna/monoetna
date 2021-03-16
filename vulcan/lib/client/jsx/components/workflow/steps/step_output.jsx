import React, {useContext} from 'react';
import Plot from 'react-plotly.js';

import ConsignmentTable from 'etna-js/plots/components/consignment/consignment_table';
import Consignment from 'etna-js/plots/models/consignment';
import Link from 'etna-js/components/link';

import {VulcanContext} from '../../../contexts/vulcan';
import {OUTPUT_COMPONENT} from '../../../models/steps';

import StepName from './step_name';

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

  // We need to extract the data from the input source.
  let rawInputData = uiStepInputDataRaw({step, pathIndex, status});

  let Component;
  switch (step.run.split('/')[1].replace('.cwl', '')) {
    case OUTPUT_COMPONENT.PLOTLY:
      // Plotly.js data payload should be in format of:
      // {
      //   data: <JSON>,
      //   layout: <JSON>
      // }
      if (!rawInputData || !rawInputData.data || !rawInputData.layout)
        return null;
      Component = (
        <Plot data={rawInputData.data} layout={rawInputData.layout}></Plot>
      );
      break;
    case OUTPUT_COMPONENT.CONSIGNMENT:
      // Not sure how to check consignment format?
      if (!rawInputData) return null;
      Component = (
        <div className='consignment-view'>
          <ConsignmentTable
            consignment={new Consignment(rawInputData)}
          ></ConsignmentTable>
        </div>
      );
      break;
    case OUTPUT_COMPONENT.RAW:
      if (!rawInputData) return null;
      Component = <div className='raw-view'>{rawInputData}</div>;
      break;
    case OUTPUT_COMPONENT.LINK:
    default:
      Component = (
        <Link link={uiStepInputDataLink({step, pathIndex, status})}>
          Download data here
        </Link>
      );
      break;
  }

  if (!Component) return null;

  return (
    <div className='step-output'>
      <StepName
        step={step}
        status={status[pathIndex][stepIndex].status}
      ></StepName>
      <div className='outputs-pane'>{Component}</div>
    </div>
  );
}
