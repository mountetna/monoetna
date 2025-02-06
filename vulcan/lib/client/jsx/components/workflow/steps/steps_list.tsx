import React, {useContext, useState, useEffect, useMemo} from 'react';
import Icon from 'etna-js/components/icon';

import {VulcanContext} from '../../../contexts/vulcan_context';
import Step from './step';
import {completedUiOutputSteps, stepOfName} from '../../../selectors/workflow_selectors';
import {useWorkspace} from '../../../contexts/workspace_context';
import { WorkspaceStep } from '../../../api_types';

export default function StepsList() {
  const [open, setOpen] = useState(false);
  const {state} = useContext(VulcanContext);
  const {status} = state;
  const {workspace} = useWorkspace();
  if (!workspace.vulcan_config) return null;

  function handleToggle() {
    setOpen(!open);
  }

  const outputs = useMemo(() => completedUiOutputSteps(workspace, status), [workspace, status]);
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
        {/* ToDo: Order based on workspace.dag */}
        {workspace ? workspace.dag.map((step, index) => {
              return <Step key={index} step={stepOfName(step, workspace.vulcan_config) as WorkspaceStep}/>;
            })
          : null}
      </div>
    </div>
  );
}
