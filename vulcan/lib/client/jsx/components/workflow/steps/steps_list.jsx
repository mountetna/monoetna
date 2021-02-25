import React, {useContext, useEffect, useState} from 'react';

import {VulcanContext} from '../../../contexts/vulcan';
import Step from './step';
import StepViewCard from './step_view_card';
import {validPath} from '../../../utils/workflow';

export default function StepsList() {
  const {workflow, pathIndex} = useContext(VulcanContext);

  const [cardIndex, setCardIndex] = useState(null);
  const [cardStep, setCardStep] = useState(null);

  function handleOnClick(index) {
    setCardIndex(cardIndex === index ? null : index);
  }

  useEffect(() => {
    if (null === cardIndex) {
      setCardStep(null);
    } else {
      setCardStep(workflow.steps[pathIndex][cardIndex]);
    }
  }, [cardIndex]);

  return (
    <div className='steps-list'>
      <div className='title'>Steps</div>
      <div className='steps-list-wrapper'>
        {validPath({workflow, pathIndex})
          ? workflow.steps[pathIndex].map((step, index) => {
              return (
                <Step
                  key={index}
                  step={step}
                  index={index}
                  onClick={handleOnClick}
                ></Step>
              );
            })
          : null}
      </div>
      <div>
        {cardStep ? (
          <StepViewCard step={cardStep} stepIndex={cardIndex}></StepViewCard>
        ) : null}
      </div>
    </div>
  );
}
