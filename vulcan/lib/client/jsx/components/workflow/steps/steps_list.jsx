import React, {useContext, useEffect, useState} from 'react';

import {VulcanContext} from '../../../contexts/vulcan';
import Step from './step';

export default function StepsList() {
  const {workflow, pathIndex, stepIndex} = useContext(VulcanContext);

  const [activeIndex, setActiveIndex] = useState(false);

  useEffect(() => {
    setActiveIndex(stepIndex);
  }, []);

  useEffect(() => {
    setActiveIndex(stepIndex);
  }, [stepIndex]);

  return (
    <ul className='steps-list'>
      {workflow.steps && workflow.steps[pathIndex]
        ? workflow.steps[pathIndex].map((step, index) => {
            return (
              <Step
                key={index}
                step={step}
                index={index}
                active={activeIndex === index}
              ></Step>
            );
          })
        : null}
    </ul>
  );
}
