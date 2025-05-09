// Framework libraries.
import React, {useContext, useEffect} from 'react';
import 'regenerator-runtime/runtime';
import {VulcanContext} from '../contexts/vulcan_context';
import {useFeatureFlag} from "etna-js/hooks/useFeatureFlag";

import WorkspaceInitializer from './workspace/workspace_initializer';
import Typography from '@material-ui/core/Typography';

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

  if (!useFeatureFlag('vulcan')) {
    return <div style={{padding: '10px'}}>
      <Typography variant='h3'>
        Vulcan is under development
      </Typography>
      <Typography>
        Access is restricted to a small group of testers at this stage, and it appears that you are not one of these testers.
      </Typography>
      <Typography>
        Contact the engineering team if you believe you are seeing this message in error.
      </Typography>
    </div>
  }

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
