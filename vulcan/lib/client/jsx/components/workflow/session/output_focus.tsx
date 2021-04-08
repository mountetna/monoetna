import React, {useEffect, useContext, useMemo} from 'react';
import {VulcanContext} from '../../../contexts/vulcan_context';
import {completedSteps} from "../../../selectors/workflow_selectors";

// TODO: Make this a useOutputFocus function instead.
export default function OutputFocus() {
  const {state} = useContext(VulcanContext);
  const {status, workflow} = state;
  if (workflow == null) return null;

  const completed = completedSteps(workflow, status);

  useEffect(() => {
    if (completed.length > 0) {
      const outputNodes = document.querySelectorAll('.step-output');
      let lastOutput;

      if (outputNodes.length > 0)
        lastOutput = outputNodes[outputNodes.length - 1];

      if (lastOutput) {
        lastOutput.scrollIntoView({behavior: 'smooth'});
      }
    }
  }, [completed]);

  // This component does not render anything.
  return null;
}
