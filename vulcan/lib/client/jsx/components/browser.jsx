// Framework libraries.
import React from 'react';
import 'regenerator-runtime/runtime';

import WorkflowManager from './workflow/workflow_manager';

export default function Browser({workflowName}) {
  return (
    <main className='vulcan-browser browser'>
      <WorkflowManager workflowName={workflowName}/>
    </main>
  );
}
