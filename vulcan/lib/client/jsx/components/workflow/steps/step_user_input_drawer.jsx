import React, {useContext} from 'react';

import {VulcanContext} from '../../../contexts/vulcan';
import UserInput from '../user_interactions/inputs/user_input';

import {uiStepType, uiStepOptions} from '../../../utils/workflow';

export default function StepUserInputDrawer({step, handleInputChange}) {
  let {workflow, session, status, pathIndex} = useContext(VulcanContext);

  // We need to unpack the grouped steps and add docs
  let stepInputs = step.in.reduce((result, input) => {
    let originalStepName = input.source[0];
    let inputName = input.source.join('/');

    // We need to fetch the original step, to see if the options
    //   data is available.
    let originalStep = workflow.steps[pathIndex].find(
      (s) => originalStepName === s.name
    );

    result.push({
      type: uiStepType(step),
      label: input.label || inputName,
      default: session.inputs[inputName] || null,
      options: uiStepOptions({step: originalStep, pathIndex, status}),
      name: inputName,
      doc: input.doc
    });
    return result;
  }, []);

  return (
    <React.Fragment>
      {stepInputs.map((input, index) => {
        return (
          <UserInput
            input={input}
            hideLabel={false}
            onChange={handleInputChange}
            key={index}
          ></UserInput>
        );
      })}
    </React.Fragment>
  );
}
