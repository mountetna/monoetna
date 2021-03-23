import React, {useContext, useState, useEffect} from 'react';

import {VulcanContext} from '../../../contexts/vulcan';
import Step from './step';
import {validPath} from '../../../utils/workflow';
import {completedUiOutputsSelector} from '../../../selectors/workflow';

export default function StepsList() {
  const [open, setOpen] = useState(true);
  const context = useContext(VulcanContext);
  let {workflow, pathIndex, status} = context;

  function handleToggle() {
    setOpen(!open);
  }

  useEffect(() => {
    // Automatically close the drawer when there are
    //   output steps.
    let outputs = completedUiOutputsSelector(context);
    if (outputs.length > 0) setOpen(false);
  }, [status]);

  return (
    <div className={`steps-list toggle-control ${open ? 'open' : 'closed'}`}>
      <div className='steps-list-header' onClick={handleToggle}>
        <div className='steps-list-toggle arrow-toggle toggle-left-right'></div>
        <div className='title'>Steps</div>
      </div>
      <div className='steps-list-positioner'>
        <div className='steps-list-wrapper'>
          {validPath({workflow, pathIndex})
            ? workflow.steps[pathIndex].map((step, index) => {
                return <Step key={index} step={step} index={index}></Step>;
              })
            : null}
        </div>
      </div>
    </div>
  );
}
