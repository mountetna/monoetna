import React, {useContext, useMemo} from 'react';

import {VulcanContext} from '../../../contexts/vulcan_context';

import StepOutput from '../steps/step_output';
import OutputFocus from './output_focus';
import {completedSteps, uiOutputOfStep} from "../../../selectors/workflow_selectors";

export default function OutputFeed() {
  // Shows stream of Output, Plots, etc.,
  //   as the session object updates.
  const {state} = useContext(VulcanContext);
  const {status, workflow} = state;
  if (!workflow) return null;

  let outputs = useMemo(
      () => completedSteps(workflow, status).filter(({step}) => !!uiOutputOfStep(step)),
      [workflow, status],
  );

  return (
    <div className='session-output-feed'>
      {outputs.map((s, index) => (
        <StepOutput key={index} step={s.step}/>
      ))}
      <OutputFocus />
    </div>
  );
}
