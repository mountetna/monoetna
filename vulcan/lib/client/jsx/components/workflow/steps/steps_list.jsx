import React, {useContext, useState, useEffect} from 'react';

import Icon from 'etna-js/components/icon';
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
    //   output steps. Open drawer if the output
    //   steps are removed...
    let outputs = completedUiOutputsSelector(context);
    setOpen(outputs.length === 0);
  }, [status]);

  return (
    <div className={`steps-list toggle-control ${open ? 'open' : 'closed'}`}>
      <div className='steps-list-header' onClick={handleToggle}>
        <Icon icon={ open ? 'angle-right' : 'angle-left' } className='toggle'/>
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
