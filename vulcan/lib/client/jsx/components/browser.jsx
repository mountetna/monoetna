// Framework libraries.
import React, {useState, useEffect} from 'react';
import 'regenerator-runtime/runtime';

import {getProjects} from 'etna-js/api/janus-api';
import {projectNameFull} from 'etna-js/utils/janus';

import {ArchimedesProvider} from '../contexts/archimedes';
import WorkflowManager from './workflow/manager';

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
      <header className='header'>
        {projectNameFull(projects, CONFIG.project_name) || CONFIG.project_name}
      </header>
      <section>
        <ArchimedesProvider>
          <WorkflowManager></WorkflowManager>
        </ArchimedesProvider>
      </section>
    </main>
  );
}
