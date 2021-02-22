// Framework libraries.
import React from 'react';
import 'regenerator-runtime/runtime';

import {VulcanProvider} from '../contexts/vulcan';
import WorkflowManager from './workflow/workflow_manager';

export default function Browser() {
  return (
    <main className='vulcan-browser browser'>
      <section>
        <VulcanProvider>
          <WorkflowManager></WorkflowManager>
        </VulcanProvider>
      </section>
    </main>
  );
}
