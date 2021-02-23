import React, {useContext, useEffect, useState} from 'react';

import {VulcanContext} from '../../../contexts/vulcan';
import Step from './step';
import {validPath} from '../../../utils/workflow';

export default function StepsList() {
  const {workflow, pathIndex} = useContext(VulcanContext);

  return (
    <div className='steps-list'>
      <div className='title'>Steps</div>
      {validPath({workflow, pathIndex})
        ? workflow.steps[pathIndex].map((step, index) => {
            return <Step key={index} step={step} index={index}></Step>;
          })
        : null}
    </div>
  );
}
