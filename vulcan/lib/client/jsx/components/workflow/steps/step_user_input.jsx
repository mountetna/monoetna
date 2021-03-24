import React, {useContext} from 'react';

import {VulcanContext} from '../../../contexts/vulcan';
import UserInput from '../user_interactions/inputs/user_input';

import {
  uiStepInputNames,
  uiStepType,
  uiStepOptions
} from '../../../utils/workflow';

export default function StepUserInput({step, handleInputChange}) {
  let {session, status, pathIndex} = useContext(VulcanContext);

  // We need to map the user input step's output to
  //   a set of input items.
  let inputNames = uiStepInputNames(step);
  let stepInputs = inputNames.reduce((result, inputName) => {
    result.push({
      type: uiStepType(step),
      label: step.label || step.name,
      default: session.inputs[inputName] || null,
      options: uiStepOptions({step, pathIndex, status}),
      name: inputName
    });
    return result;
  }, []);

  return (
    <React.Fragment>
      {stepInputs.map((input, index) => {
        return (
          <UserInput
            input={input}
            hideLabel={true}
            onChange={handleInputChange}
            key={index}
          ></UserInput>
        );
      })}
    </React.Fragment>
  );
}
