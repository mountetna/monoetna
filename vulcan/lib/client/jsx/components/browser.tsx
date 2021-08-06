// Framework libraries.
import React from 'react';
import 'regenerator-runtime/runtime';

import WorkflowManager from './workflow/workflow_manager';

interface Props {
  workflowName: string,
  projectName: string,
}

export default function Browser({workflowName, projectName}: Props) {
  return (
    <main className='vulcan-browser browser'>
      <WorkflowManager workflowName={workflowName} projectName={projectName}/>
    </main>
  );
}
