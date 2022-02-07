// Framework libraries.
import React from 'react';
import 'regenerator-runtime/runtime';

import WorkflowManager from './workflow/workflow_manager';

interface Props {
  workflowName: string,
  projectName: string,
}

export default function Browser({workflow_name, project_name, figure_id}: Props) {
  return (
    <main className='vulcan-browser browser'>
      <WorkflowManager
        workflowName={workflow_name}
        figureId={figure_id}
        projectName={project_name}/>
    </main>
  );
}
