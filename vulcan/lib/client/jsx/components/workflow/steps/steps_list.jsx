import React, {useContext, useEffect, useState} from 'react';

import {VulcanContext} from '../../../contexts/vulcan';
import Step from './step';
import {validPath} from '../../../selectors/workflow_selector';

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
    <div className='steps-list'>
      <div className='title'>Steps</div>
      {validPath({workflow, pathIndex})
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
    </div>
  );
}
