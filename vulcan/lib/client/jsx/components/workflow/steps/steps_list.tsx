import React, {useContext, useState, useEffect, useMemo} from 'react';
import Icon from 'etna-js/components/icon';

import {VulcanContext} from '../../../contexts/vulcan_context';
import Step from './step';
import {completedUiOutputSteps} from "../../../selectors/workflow_selectors";
import {useWorkflow} from "../../../contexts/workflow_context";

export default function StepsList() {
  const [open, setOpen] = useState(false);
  const {state} = useContext(VulcanContext);
  const {status} = state;
  const {workflow} = useWorkflow();

  function handleToggle() {
    setOpen(!open);
  }

  const outputs = useMemo(() => completedUiOutputSteps(workflow, status), [workflow, status]);
  const hasCompletedOutputs = outputs.length === 0;

  useEffect(() => {
    setOpen(hasCompletedOutputs || !!state.pollingState);
  }, [state.pollingState, hasCompletedOutputs]);

  return (
    <div className={`steps-list toggle-control ${open ? 'open' : 'closed'}`}>
      <div className='steps-list-header' onClick={handleToggle}>
        <Icon icon={ open ? 'angle-right' : 'angle-left' } className='toggle'/>
        <div className='title'>Steps</div>
      </div>
      <div className='steps-list-wrapper'>
        {workflow ? workflow.steps[0].map((step, index) => {
              return <Step key={index} step={step}/>;
            })
          : null}
      </div>
    </div>
  );
}
