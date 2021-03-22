import React, {useContext, useState} from 'react';

import {VulcanContext} from '../../../contexts/vulcan';
import Step from './step';
import {validPath} from '../../../utils/workflow';

export default function StepsList() {
  const [open, setOpen] = useState(false);
  const {workflow, pathIndex} = useContext(VulcanContext);

  function handleToggle() {
    setOpen(!open);
  }

  return (
    <div className={`steps-list toggle-control ${open ? 'open' : 'closed'}`}>
      <div className='steps-list-header' onClick={handleToggle}>
        <div className='steps-list-toggle arrow-toggle toggle-left-right'></div>
        <div className='title'>Steps</div>
      </div>
      <div className='steps-list-wrapper'>
        {validPath({workflow, pathIndex})
          ? workflow.steps[pathIndex].map((step, index) => {
              return <Step key={index} step={step} index={index}></Step>;
            })
          : null}
      </div>
    </div>
  );
}
