// Framework libraries.
import React, {useContext, useEffect} from 'react';
import 'regenerator-runtime/runtime';
import {VulcanContext} from '../contexts/vulcan_context';

import WorkspaceInitializer from './workflow/workspace_initializer';

interface Props {
  workflowName: string;
  projectName: string;
}

export default function Browser({
  workflow_name,
  project_name,
  workspace_id
}: {
  workflow_name: string;
  project_name: string;
  workspace_id: number;
}) {
  const {
    state: {workflows}
  } = useContext(VulcanContext);

  if (workflows.length === 0 || !project_name) return null;

  return (
    <main className='vulcan-browser browser'>
      <WorkspaceInitializer
        workflowName={workflow_name}
        workspaceId={workspace_id}
        projectName={project_name}
      />
    </main>
  );
}
