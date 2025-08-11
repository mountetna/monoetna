import React, {useContext, useState, useEffect, useMemo} from 'react';
import Icon from 'etna-js/components/icon';

import {VulcanContext} from '../../contexts/vulcan_context';
import {outputUINames, outputUIsWithInputsReady, stepOfName} from '../../selectors/workflow_selectors';
import {useWorkspace} from '../../contexts/workspace_context';
import { VulcanConfig } from '../../api_types';
import StepIconName from './drawers/step_elements/step_icon_name';
import FlatButton from 'etna-js/components/flat-button'
import Grid from '@material-ui/core/Grid';

export default function StepsList() {
  const [open, setOpen] = useState(false);
  const {state} = useContext(VulcanContext);
  const {status, workspace, refreshingProgress} = state;

  function handleToggle() {
    setOpen(!open);
  }

  const hasCompletedOutputs = useMemo(() => {
    if (!workspace) return false;
    return outputUIsWithInputsReady(workspace, status.file_contents).length > 0
  }, [workspace, status.file_contents]);

  useEffect(() => {
    setOpen(!hasCompletedOutputs || state.isSyncing);
  }, [state.isSyncing, hasCompletedOutputs]);

  if (!workspace) return null;
  const stepNamesToStepIconNames = (step, index) => (
    <div key={index} className='step'>
      <StepIconName step={stepOfName(step, workspace.vulcan_config)}/>
    </div>
  )

  return (
    <div className={`steps-list toggle-control ${open ? 'open' : 'closed'}`}>
      <div className='steps-list-header' onClick={handleToggle}>
        <Icon icon={ open ? 'angle-right' : 'angle-left' } className='toggle'/>
        {open ? <div className='title'>Progress</div> : null}
      </div>
      <div className='steps-list-subheader'>Steps</div>
      <div className='steps-list-wrapper'>
        {workspace.dag.map(stepNamesToStepIconNames)}
      </div>
      <Grid className='steps-list-subheader' container direction='row' justifyContent='space-between' alignItems='flex-end'>
        <Grid item>Output</Grid>
        <Grid item><FlatButton
          className='header-btn-name-save'
          icon={false? 'spinner fa-spin' : 'circle-notch'}
          label='Refresh'
          title={false? 'Coming Soon: Assess Output Validity' : 'Coming Soon: Assessing Output Validity'}
          disabled={true}
          onClick={() => {
          }}
        /></Grid>
      </Grid>
      <div className='steps-list-wrapper'>
        {outputUINames(workspace).map(stepNamesToStepIconNames)}
      </div>
    </div>
  );
};
