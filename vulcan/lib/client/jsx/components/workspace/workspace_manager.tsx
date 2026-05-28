import React, {
  useCallback,
  useEffect,
  useState,
  useContext,
  useMemo
} from 'react';

import {VulcanContext} from '../../contexts/vulcan_context';
import {WorkspaceProvider} from '../../contexts/workspace_context';
import InputFeed from './input_feed';
import OutputFeed from './output_feed';
import StepsList from './steps_list';

import WorkspaceBreadcrumb from './workspace_breadcrumb';
import WorkspaceManagerControls from './workspace_manager_controls';

const WorkspaceManager = ({workspaceId}) => {
  const [saving, setSaving] = useState(false);
  const { state } = useContext(VulcanContext);
  
  const {projectName, workflows, workspaces} = state;

  const workspace = workspaces.find( s => s.workspace_id == workspaceId );

  const workflow = workflows.find( f => f.id == workspace?.workflow_id);

  if (!workflow || !workspace) return <div>No such workspace {workspaceId}</div>;

  return (
    <WorkspaceProvider projectName={ projectName} workspace={workspace} workflow={workflow}>
      <div className='session-manager'>
        <div className='session-header'>
          <WorkspaceBreadcrumb/>
          <WorkspaceManagerControls/>
        </div>
        <div className='session-feed-container'>
          <OutputFeed/>
          <StepsList/>
        </div>
      </div>
    </WorkspaceProvider>
  );
}

export default WorkspaceManager;
