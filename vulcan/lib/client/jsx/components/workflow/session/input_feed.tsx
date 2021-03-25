import React, {useContext, useEffect} from 'react';

import {VulcanContext} from '../../../contexts/vulcan_context';

import PrimaryInputs from './primary_inputs';
import StepError from '../steps/step_error';
import {
  completedSteps, erroredSteps, groupUiSteps,
  isPendingUiQuery,
  missingUiQueryOutputs,
  pendingSteps,
  uiQueryOfStep
} from "../../../selectors/workflow_selectors";
import {setInputs} from "../../../actions/vulcan";
import StepUserInputWrapper from "../steps/step_user_input_wrapper";

export default function InputFeed() {
  // Shows stream of Inputs,
  //   as the session object updates.
  const {state, dispatch} = useContext(VulcanContext);
  const {workflow, session, status, data, inputs} = state;

  if (!workflow) return null;

  let completed = completedSteps(workflow, status).filter(({step}) => !!uiQueryOfStep(step));
  let nextUiSteps = pendingSteps(workflow, status).filter(({step}) => isPendingUiQuery(step, status, data, session));

  // We inject a `null` input into the session,
  //   to indicate that we're waiting for a user input value
  //   to return to the server.
  useEffect(() => {
    let newInputs = inputs;

    nextUiSteps.forEach((nextInputStep) => {
      let missingInputs = missingUiQueryOutputs(nextInputStep.step, inputs);

      if (Object.keys(missingInputs).length > 0) {
        // Make sure to copy over the current inputs, otherwise
        //   they'll get wiped out in the reducer.
        newInputs = {
          ...newInputs,
          ...missingInputs,
        };
      }
    });

    if (newInputs !== inputs) dispatch(setInputs(newInputs));
  }, [nextUiSteps, inputs]);

  const uiSteps = groupUiSteps(completed.concat(nextUiSteps));

  let errorSteps = erroredSteps(workflow, status);

  return (
    <div className='session-input-feed'>
      <PrimaryInputs/>
      {uiSteps.map((s, index) => (
        <StepUserInputWrapper
          key={index}
          step={s.step}
        />
      ))}
      {errorSteps.map((s, index) => (
        <StepError key={index} step={s.step}/>
      ))}
    </div>
  );
}
