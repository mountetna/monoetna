import React, {useState, useContext, useEffect} from 'react';

import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
import {showMessages} from 'etna-js/actions/message_actions';

import Dropdown from 'etna-js/components/inputs/dropdown';

import {VulcanContext} from '../../contexts/vulcan';
import {getWorkflow, getWorkflows} from '../../api/vulcan';
import StepsList from './steps/steps_list';
import CurrentStep from './steps/current_step';

export default function Input() {
  const invoke = useActionInvoker();
  const {
    workflow,
    pathIndex,
    stepIndex,
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
    <section className='input-section'>
      <div className='steps-list'>
        <StepsList></StepsList>
      </div>
      <div className='steps-current'>
        <CurrentStep></CurrentStep>
      </div>
    </section>
  );
}
