import React, {useContext, useState, useEffect} from 'react';
import Icon from 'etna-js/components/icon';

import {VulcanContext} from '../../../contexts/vulcan_context';
import Step from './step';
import {completedUiOutputSteps} from "../../../selectors/workflow_selectors";

export default function StepsList() {
  const [open, setOpen] = useState(false);
  const {state} = useContext(VulcanContext);
  const {workflow, status} = state;

  function handleToggle() {
    setOpen(!open);
  }

  useEffect(() => {
    // Automatically close the drawer when there are
    //   output steps. Open drawer if the output
    //   steps are removed...
    if (!workflow) return;
    let outputs = completedUiOutputSteps(workflow, status);
    setOpen(outputs.length === 0);
  }, [status]);

  return (
    <div className={`steps-list toggle-control ${open ? 'open' : 'closed'}`}>
      <div className='steps-list-header' onClick={handleToggle}>
        <Icon icon={ open ? 'angle-right' : 'angle-left' } className='toggle'/>
        <div className='title'>Steps</div>
      </div>
      <div className='steps-list-wrapper'>
        {workflow ? workflow.steps[0].map((step, index) => {
              return <Step key={index} step={step}></Step>;
            })
          : null}
      </div>
    </div>
  );
}
