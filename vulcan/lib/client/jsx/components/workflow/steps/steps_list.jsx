import React, {useContext} from 'react';

import {VulcanContext} from '../../../contexts/vulcan';
import Step from './step';

export default function StepsList() {
  const {workflow, pathIndex, setStepIndex} = useContext(VulcanContext);

  function handleOnClick(index) {
    setStepIndex(index);
  }

  return (
    <ol className='steps-list'>
      {workflow.steps && workflow.steps[pathIndex]
        ? workflow.steps[pathIndex].forEach((step, index) => {
            return (
              <Step
                key={index}
                step={step}
                index={index}
                onClick={(e) => handleOnClick(index)}
              ></Step>
            );
          })
        : null}
    </ol>
  );
}
