// Framework libraries.
import React, {useState, useCallback, useEffect, useMemo} from 'react';
import 'regenerator-runtime/runtime';
import useAsyncWork from 'etna-js/hooks/useAsyncWork';
import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';

import {getProjects} from 'etna-js/api/janus-api';
import Dropdown from 'etna-js/components/inputs/dropdown';

// Module imports.
import {setLocation} from 'etna-js/actions/location_actions';
import {showMessages} from 'etna-js/actions/message_actions';

import {projectNameFull} from 'etna-js/utils/janus';
import {getWorkflow, getWorkflows} from '../../api/archimedes_api';

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
  const invoke = useActionInvoker();
  const [workflows, setWorkflows] = useState(null);
  const [selectedWorkflowName, setSelectedWorkflowName] = useState(null);
  const [currentWorkflow, setCurrentWorkflow] = useState(null);
  const [projects, setProjects] = useState(null);

  useEffect(() => {
    getProjects().then((projects) => {
      setProjects(projects);
    });
    getWorkflows()
      .then((allWorkflows) => {
        console.log('all workflows');
        console.log(allWorkflows);
        setWorkflows(allWorkflows);
      })
      .catch((e) => {
        invoke(showMessages(e));
      });
  }, []);

  useEffect(() => {
    if (selectedWorkflowName) {
      getWorkflow(selectedWorkflowName)
        .then((workflowDetails) => {
          console.log(workflowDetails);
          setCurrentWorkflow(workflowDetails);
        })
        .catch((e) => {
          invoke(showMessages(e));
        });
    }
  }, [selectedWorkflowName]);

  return (
    <div className='vulcan-browser'>
      <div>
        Select a workflow for{' '}
        {projectNameFull(projects, CONFIG.project_name) || CONFIG.project_name}
      </div>
      <div>
        <Dropdown
          list={workflows ? Object.keys(workflows) : []}
          default_text='Select a workflow'
          onSelect={(e) => setSelectedWorkflowName(e)}
        ></Dropdown>
      </div>
    </div>
  );
}

function useWorkflowActions(setCurrentWorkflow) {
  const invoke = useActionInvoker();
  const {revision, model_name, template, record_name} = browserState;

  return {
    cancelEdits,
    approveEdits
  };

  function cancelEdits() {
    setMode('browse');
    invoke(discardRevision(record_name, model_name));
  }

  function postEdits() {
    setMode('submit');

    invoke(
      sendRevisions(
        model_name,
        template,
        {[record_name]: revision},
        () => setMode('browse'),
        () => setMode('edit')
      )
    );
  }

  function approveEdits() {
    if (Object.keys(revision).length > 0) postEdits();
    else cancelEdits();
  }
}
