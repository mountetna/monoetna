import React, {useContext, useState} from 'react';

import {VulcanContext} from '../../../contexts/vulcan_context';
import Step from './step';

export default function StepsList() {
  const [open, setOpen] = useState(false);
  const {state} = useContext(VulcanContext);
  const {workflow} = state;

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
        { workflow
          ? workflow.steps[0].map((step) => {
              return <Step key={step.name} step={step}/>;
            })
          : null}
      </div>
    </div>
  );
}
