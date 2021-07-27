import React, {useContext, useEffect} from 'react';

import {VulcanContext} from '../../../contexts/vulcan_context';

import PrimaryInputs from './primary_inputs';
import StepError from '../steps/step_error';
import {
  completedSteps, erroredSteps, groupUiSteps,
  isPendingUiQuery,
  pendingSteps,
  uiQueryOfStep
} from "../../../selectors/workflow_selectors";
import StepUserInputWrapper from "../steps/step_user_input_wrapper";

export default function InputFeed() {
  // Shows stream of Inputs,
  //   as the session object updates.
  const {state} = useContext(VulcanContext);
  const {workflow, session, status, data} = state;

  if (!workflow) return null;

  let completed = completedSteps(workflow, status).filter(step => !!uiQueryOfStep(step));
  let nextUiSteps = pendingSteps(workflow, status).filter(step => isPendingUiQuery(step, status, data, session));
  const groupedSteps = groupUiSteps(completed.concat(nextUiSteps));

  let errorSteps = erroredSteps(workflow, status);

  return (
    <div className='session-input-feed'>
      <PrimaryInputs/>
      {groupedSteps.map((s, index) => (
        <StepUserInputWrapper
          key={index}
          group={s}
        />
      ))}
      {errorSteps.map((s, index) => (
        <StepError key={index} step={s.step}/>
      ))}
    </div>
  );
}
