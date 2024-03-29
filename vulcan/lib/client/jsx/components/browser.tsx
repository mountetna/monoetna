// Framework libraries.
import React, {useContext} from 'react';
import 'regenerator-runtime/runtime';
import {VulcanContext} from '../contexts/vulcan_context';

import WorkflowManager from './workflow/workflow_manager';

interface Props {
  workflowName: string;
  projectName: string;
}

export default function Browser({
  workflow_name,
  project_name,
  figure_id
}: {
  workflow_name: string;
  project_name: string;
  figure_id: number;
}) {
  const {
    state: {workflows}
  } = useContext(VulcanContext);

  if (workflows.length === 0 || !project_name) return null;

  return (
    <main className='vulcan-browser browser'>
      <WorkflowManager
        workflowName={workflow_name}
        figureId={figure_id}
        projectName={project_name}
      />
    </main>
  );
}
