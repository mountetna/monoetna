import React, {useContext, useEffect, useState} from 'react';

import {VulcanContext} from '../../../contexts/vulcan';
import {STATUS} from '../../../models/steps';
import {validStep} from '../../../selectors/workflow_selector';

import StepComplete from './step_complete';
import StepError from './step_error';
import StepPending from './step_pending';

export default function ActiveStep() {
  const {workflow, pathIndex, stepIndex, status} = useContext(VulcanContext);

  if (!validStep({workflow, pathIndex, stepIndex})) return null;

  let activeStepStatus = status[pathIndex][stepIndex];

  let Component;

  switch (activeStepStatus.status) {
    case STATUS.COMPLETE:
      Component = StepComplete;
      break;
    case STATUS.ERROR:
      Component = StepError;
      break;
    case STATUS.PENDING:
      Component = StepPending;
      break;
    default:
      break;
  }

  if (!Component) return null;

  return <Component></Component>;
}
