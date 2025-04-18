import React, {useContext, useMemo} from 'react';

import {VulcanContext} from '../../../contexts/vulcan_context';

import StepOutput from '../steps/step_output';
import {completedUiOutputSteps} from '../../../selectors/workflow_selectors';
import {useWorkflow} from '../../../contexts/workflow_context';

export default function OutputFeed() {
  // Shows stream of Output, Plots, etc.,
  //   as the session object updates.
  const {state} = useContext(VulcanContext);
  const {workflow} = useWorkflow();
  const {status} = state;

  let outputs = useMemo(
      () => completedUiOutputSteps(workflow, status),
      [workflow, status],
  );

  return (
    <div className="session-output-feed">
      {outputs.map((s, index) => (
        <StepOutput key={index} step={s}/>
      ))}
    </div>
  );
}
