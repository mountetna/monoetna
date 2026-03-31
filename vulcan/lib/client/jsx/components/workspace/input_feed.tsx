import React, {useContext, useEffect, useMemo} from 'react';

import {VulcanContext} from '../../contexts/vulcan_context';

import ParamInputs from './drawers/param_user_inputs';
import StepError from './drawers/step_error';
import {
  completedStepNames,
  erroredStepNames,
  groupUiSteps,
  pendingUIInputStepReady,
  pendingStepNames,
  uiComponentOfStep,
  stepNamesOfStatus,
  pickToArray,
  stepOfName
} from '../../selectors/workflow_selectors';
import GroupedStepUI from './drawers/step_user_input';
import { LoadingIconWithText } from '../dashboard/loading_icon';
import { makeStyles } from '@material-ui/core/styles';

export const useInputFeedStyles = makeStyles((theme) => ({
  card: {
    borderRadius: 0,
    border: '1px solid #eee'
  },
  error: {
    border: '1px solid red'
  },
  label: {
    paddingLeft: '5px'
  },
  header: {
    cursor: 'pointer',
    padding: '10px'
  }
}));

export default function InputFeed() {
  // Shows stream of Inputs,
  //   as the session object updates.
  const {state} = useContext(VulcanContext);
  const {workflow, workspace, status, update_files} = state;

  if (!workflow.name || !workspace || !workspace.vulcan_config) return null;

  let completed = completedStepNames(workspace, status).filter(
    (step) => !!uiComponentOfStep(step, workspace.vulcan_config)
  );
  let nextUiSteps = stepNamesOfStatus(['pending', 'upcoming'], workspace, status).filter((step) =>
    pendingUIInputStepReady(step, status, workspace)
  );
  const groupedSteps = groupUiSteps(completed.concat(nextUiSteps), workspace);

  const stepInputs = useMemo(() => {
    return update_files ?
      <LoadingIconWithText text='Refreshing Files'/>:
      groupedSteps.map((s, index) => (
        <GroupedStepUI key={index} group={s} />
      ));
  }, [update_files, status.file_contents])

  let errorSteps = erroredStepNames(workspace, status).map((step) => stepOfName(step, workspace.vulcan_config));

  return (
    <div className='session-input-feed'>
      <ParamInputs />
      {stepInputs}
      {errorSteps.map((s, index) => (
        <StepError key={index} step={s} />
      ))}
    </div>
  );
}
