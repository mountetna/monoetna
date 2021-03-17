import React, {useContext} from 'react';
import Plot from 'react-plotly.js';

import ConsignmentTable from 'etna-js/plots/components/consignment/consignment_table';
import Consignment from 'etna-js/plots/models/consignment';
import Link from 'etna-js/components/link';

import {VulcanContext} from '../../../contexts/vulcan';
import {OUTPUT_COMPONENT} from '../../../models/steps';

import StepName from './step_name';
import UserOutput from '../user_interactions/user_output';

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

  return (
    <div className='step-output'>
      <StepName
        step={step}
        status={status[pathIndex][stepIndex].status}
      ></StepName>
      <div className='outputs-pane'>
        <UserOutput step={step}></UserOutput>
      </div>
    </div>
  );
}
