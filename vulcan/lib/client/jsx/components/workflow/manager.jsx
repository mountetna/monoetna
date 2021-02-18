import React, {useState, useContext, useEffect} from 'react';

import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
import {showMessages} from 'etna-js/actions/message_actions';

import Dropdown from 'etna-js/components/inputs/dropdown';

import {VulcanContext} from '../../contexts/vulcan';
import {getWorkflow, getWorkflows} from '../../api/vulcan';
import Output from './output';

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

export default function Manager() {
  const invoke = useActionInvoker();
  const {
    workflows,
    workflow,
    pathIndex,
    stepIndex,
    setWorkflows,
    setWorkflow,
    setPathIndex,
    setStepIndex
  } = useContext(VulcanContext);

  const [selectedWorkflowName, setSelectedWorkflowName] = useState(null);

  useEffect(() => {
    getWorkflows()
      .then((allWorkflows) => {
        setWorkflows(allWorkflows);
      })
      .catch((e) => {
        invoke(showMessages(e));
      });
  }, []);

  useEffect(() => {
    if (selectedWorkflowName) {
      getWorkflow(selectedWorkflowName, true)
        .then((workflowDetails) => {
          console.log(workflowDetails);
          setWorkflow(workflowDetails);
        })
        .catch((e) => {
          invoke(showMessages(e));
        });
    }
  }, [selectedWorkflowName]);

  function handleOnSelect(e) {
    setSelectedWorkflowName(Object.keys(workflows)[e]);
  }

  return (
    <div className='workflow-manager'>
      <Dropdown
        list={workflows ? Object.keys(workflows) : []}
        default_text='Select a workflow'
        onSelect={handleOnSelect}
      ></Dropdown>
      <div>You selected to start workflow: {selectedWorkflowName}.</div>
      <div>
        Steps:
        <ol>
          {/* workflow.steps is actually an Array of Arrays.
           This example assumes there is only one nested Array, because
           there is only one possible path. For more complicated
           workflows, we'll have to figure out how to render those. */}
          {workflow && workflow.steps
            ? workflow.steps[0].map((step, index) => {
                return <li key={index}>{step.name}</li>;
              })
            : "You'll see a list of steps here once you select a workflow."}
        </ol>
      </div>
      <Output></Output>
    </div>
  );
}
