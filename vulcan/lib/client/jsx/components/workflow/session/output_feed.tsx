import React, {useContext, useMemo} from 'react';

import {VulcanContext} from '../../../contexts/vulcan_context';

import StepOutput from '../steps/step_output';
import { outputUIsWithInputsReady } from '../../../selectors/workflow_selectors';
import {useWorkspace} from '../../../contexts/workspace_context';

export default function OutputFeed() {
  // Shows stream of Output, Plots, etc.,
  //   as the session object updates.
  const {state} = useContext(VulcanContext);
  const {workspace} = useWorkspace();
  const {status} = state;

  let outputs = useMemo(
      () => outputUIsWithInputsReady(workspace, status),
      [workspace, status],
  );

  return (
    <div className="session-output-feed">
      {outputs.map((s, index) => (
        <StepOutput key={index} step={s}/>
      ))}
    </div>
  );
}
