// Framework libraries.
import React, {useState, useEffect} from 'react';
import 'regenerator-runtime/runtime';

import {getProjects} from 'etna-js/api/janus-api';
import {projectNameFull} from 'etna-js/utils/janus';

import {VulcanProvider} from '../contexts/vulcan';
import WorkflowManager from './workflow/workflow_manager';

const loadingDiv = (
  <div className='browser'>
    <div id='loader-container'>
      <div className='loader'>Loading...</div>
    </div>
  </div>
);

const errorDiv = (
  <div className='browser'>
    <div id='loader-container'>
      <div className='loader'>Failed to load.</div>
    </div>
  </div>
);

export default function Browser() {
  const [projects, setProjects] = useState(null);

  useEffect(() => {
    getProjects().then((projects) => {
      setProjects(projects);
    });
  }, []);

  return (
    <main className='vulcan-browser browser'>
      <VulcanProvider>
        <WorkflowManager></WorkflowManager>
      </VulcanProvider>
    </main>
  );
}
