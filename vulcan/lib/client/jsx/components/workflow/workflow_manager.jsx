import React, {useContext, useEffect} from 'react';

import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
import {showMessages} from 'etna-js/actions/message_actions';

import {VulcanContext} from '../../contexts/vulcan';
import {getWorkflows} from '../../api/vulcan';
import Output from './output';
import Input from './input';

import SessionManager from './session/session_manager';

import StepsList from './steps/steps_list';

// Hardcode for now, since only one workflow
const WORKFLOW_NAME = 'umap.cwl';

export default function WorkflowManager() {
  const invoke = useActionInvoker();
  const {setWorkflow, setPathIndex} = useContext(VulcanContext);

  useEffect(() => {
    getWorkflows()
      .then((response) => {
        let currentWorkflow = response.workflows.find(
          (w) => WORKFLOW_NAME === w.name
        );
        setWorkflow(currentWorkflow);

        // longest step chain == default path for now?
        setPathIndex(
          currentWorkflow.steps
            .map((a) => a.length)
            .indexOf(Math.max(...currentWorkflow.steps.map((a) => a.length)))
        );
      })
      .catch((e) => {
        console.error(e);
        invoke(showMessages([e]));
      });
  }, []);

  return (
    <div className='workflow-manager'>
      <div className='workflow-header'>
        You are working on a {WORKFLOW_NAME.toUpperCase()} workflow:
      </div>
      <div className='step-wrapper'>
        <div className='step-nav-wrapper'>
          <StepsList></StepsList>
        </div>
        <div className='step-main-pane-wrapper'>
          <SessionManager></SessionManager>
        </div>
      </div>
    </div>
  );
}
