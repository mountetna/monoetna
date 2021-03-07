import React, {useEffect, useContext} from 'react';
import {VulcanContext} from '../../../contexts/vulcan';
import {completedUiOutputsSelector} from '../../../selectors/workflow';

export default function OutputFocus() {
  const context = useContext(VulcanContext);

  let outputs = completedUiOutputsSelector(context);

  useEffect(() => {
    if (outputs.length > 0) {
      const outputNodes = document.querySelectorAll('.step-output');
      let lastOutput;

      if (outputNodes.length > 0)
        lastOutput = outputNodes[outputNodes.length - 1];

      if (lastOutput) {
        lastOutput.scrollIntoView({behavior: 'smooth'});
      }
    }
  }, [context]);

  // This component does not render anything.
  return null;
}
