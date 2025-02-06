import React, {useContext, useEffect} from 'react';

import {VulcanContext} from '../../../contexts/vulcan_context';

import PrimaryInputs from './primary_inputs';
import StepError from '../steps/step_error';
import {
  completedStepNames,
  erroredStepNames,
  groupUiSteps,
  pendingUIInputStepReady,
  pendingStepNames,
  uiComponentOfStep,
  stepNamesOfStatus,
  pickToArray
} from '../../../selectors/workflow_selectors';
import StepUserInputWrapper from '../steps/step_user_input_wrapper';

export default function InputFeed() {
  // Shows stream of Inputs,
  //   as the session object updates.
  const {state} = useContext(VulcanContext);
  const {workflow, workspace, status,} = state;

  if (!workflow.name || !workspace || !workspace.vulcan_config) return null;

  let completed = completedStepNames(workspace, status).filter(
    (step) => !!uiComponentOfStep(step, workspace.vulcan_config)
  );
  let nextUiSteps = stepNamesOfStatus(['pending', 'upcoming'], workspace, status).filter((step) =>
    pendingUIInputStepReady(step, status, workspace)
  );
  const groupedSteps = groupUiSteps(completed.concat(nextUiSteps), workspace);

  // let errorSteps = pickToArray(workspace.steps, erroredStepNames(workspace, status));

  return (
    <div className='session-input-feed'>
      <PrimaryInputs />
      {groupedSteps.map((s, index) => (
        <StepUserInputWrapper key={index} group={s} />
      ))}
      {/* {errorSteps.map((s, index) => (
        <StepError key={index} step={s} />
      ))} */}
    </div>
  );
}
