import React, {useState, useContext, useEffect} from 'react';

import {VulcanContext} from '../../contexts/vulcan';
import StepInput from './steps/step_input';
import StepError from './steps/step_error';
import StepPending from './steps/step_pending';

import {STATUS} from '../../models/steps';
import {validStep} from '../../selectors/workflow_selector';

export default function Input() {
  const {workflow, pathIndex, stepIndex} = useContext(VulcanContext);

  if (!validStep({workflow, pathIndex, stepIndex})) return null;

  let currentStep = workflow.steps[pathIndex][stepIndex];

  let Component;
  switch (currentStep.status) {
    case STATUS.PENDING:
      Component = StepPending;
      break;
    case STATUS.ERROR:
      Component = StepError;
      break;
    case STATUS.COMPLETE:
      // Show what the inputs were at this step
      Component = StepInput;
      break;
    default:
      break;
  }

  if (!Component) return null;

  return (
    <section className='step-input'>
      <Component step={currentStep}></Component>
    </section>
  );
}
