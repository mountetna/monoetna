// Framework libraries.
import React, {useContext, useEffect} from 'react';
import 'regenerator-runtime/runtime';
import {VulcanContext} from '../contexts/vulcan_context';

import WorkspaceInitializer from './workspace/workspace_initializer';

interface Props {
  workflowName: string;
  projectName: string;
}

export default function Browser({
  project_name,
  workspace_id
}: {
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
        workspaceId={workspace_id}
        projectName={project_name}
      />
    </main>
  );
}
