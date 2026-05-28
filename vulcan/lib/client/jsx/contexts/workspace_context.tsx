import React, {useContext, useMemo} from 'react';
import {VulcanContext} from './vulcan_context';
import { defaultWorkflow, defaultWorkspace, StepStatus } from '../api_types';
import { getIsRunning, getWorkspace, pullRunStatus, getConfig } from './api';

export function useWorkspace() {
  const {state} = useContext(VulcanContext);
  const workflow = state.workflow || defaultWorkflow;
  const workspace = state.workspace || defaultWorkspace;
  const workspaceId = workspace.workspace_id;
  const {stepsStatus} = state.status.steps;

  const hasPendingEdits = state.bufferedSteps.length > 0;

  const complete = useMemo(
    () => !!stepsStatus && Object.keys(stepsStatus).length > 0 && Object.values(stepsStatus).every(({step}) => step.status === 'complete'),
    [stepsStatus]
  );

  return {workflow, workspace, workspaceId, hasPendingEdits, complete};
}

export const WorkspaceContext = React.createContext({});

export const WorkspaceProvider = ({projectName, workspace, workflow, children}) => {

  const [ isRunning, setIsRunning ] = React.useState(false);
  const [ workspaceDetails, setWorkspaceDetails] = React.useState({});
  const [ config, setConfig ] = React.useState({});
  const [ runStatus, setRunStatus ] = React.useState(null);

  const [ fileContents, setFileContents ] = React.useState({});

  const loadWorkspace = async workspace => {
    const { running } = await getIsRunning(projectName, workspace.workspace_id);
    setIsRunning(running);

    const workspaceDetails = await getWorkspace(projectName, workspace.workspace_id);
    setWorkspaceDetails(workspaceDetails);

    const newConfig = await getConfig(projectName, workspace.workspace_id);
    setConfig(newConfig);

    if (workspaceDetails.last_run_id) {
      const runStatus = await pullRunStatus( projectName, workspace.workspace_id, workspaceDetails.last_run_id )
      setRunStatus(runStatus);
    }
  }

  React.useEffect( () => {
    loadWorkspace(workspace);
  }, [workspace] );

  return <WorkspaceContext.Provider value={{
    state: {
      workspace,
      workflow,
      workspaceDetails,
      dag,
      runStatus,
      file_contents: fileContents,
      isRunning
    }
  }}>
  {children}
  </WorkspaceContext.Provider>
}
