import React, {useContext, useEffect} from 'react';

import {VulcanContext} from '../../contexts/vulcan_context';

import ParamInputs from './drawers/param_user_inputs';
import StepError from './drawers/step_error';
import {
  completedStepNames,
  erroredStepNames,
  groupUiSteps,
  pendingUIInputStepReady,
  pendingStepNames,
  uiComponentOfStep,
  stepNamesOfStatus,
  pickToArray
} from '../../selectors/workflow_selectors';
import GroupedStepUI from './drawers/step_user_input';

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
      <ParamInputs />
      {groupedSteps.map((s, index) => (
        <GroupedStepUI key={index} group={s} />
      ))}
      {/* {errorSteps.map((s, index) => (
        <StepError key={index} step={s} />
      ))} */}
    </div>
  );
}
