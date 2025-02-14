import React, {useContext, useState, useEffect, useMemo} from 'react';
import Icon from 'etna-js/components/icon';

import {VulcanContext} from '../../contexts/vulcan_context';
import {outputUIsWithInputsReady, stepOfName} from '../../selectors/workflow_selectors';
import {useWorkspace} from '../../contexts/workspace_context';
import { WorkspaceStep } from '../../api_types';
import StepIconName from './drawers/step_elements/step_icon_name';

export default function StepsList() {
  const [open, setOpen] = useState(false);
  const {state} = useContext(VulcanContext);
  const {status, workspace} = state;

  function handleToggle() {
    setOpen(!open);
  }

  const hasCompletedOutputs = useMemo(() => {
    if (!workspace) return false;
    return outputUIsWithInputsReady(workspace, status.file_contents).length > 0
  }, [workspace, status.file_contents]);

  useEffect(() => {
    setOpen(!hasCompletedOutputs || !!state.pollingState);
  }, [state.pollingState, hasCompletedOutputs]);

  return (
    <div className={`steps-list toggle-control ${open ? 'open' : 'closed'}`}>
      <div className='steps-list-header' onClick={handleToggle}>
        <Icon icon={ open ? 'angle-right' : 'angle-left' } className='toggle'/>
        <div className='title'>Steps</div>
      </div>
      <div className='steps-list-wrapper'>
        {!workspace ? null :workspace.dag.map((step, index) => (
          <div key={index} className='step'>
            <StepIconName step={stepOfName(step, workspace.vulcan_config)}/>
          </div>
        ))}
      </div>
    </div>
  );
}
