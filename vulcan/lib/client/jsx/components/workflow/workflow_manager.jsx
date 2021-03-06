import React, {useContext, useEffect} from 'react';

import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
import {showMessages} from 'etna-js/actions/message_actions';

import {VulcanContext} from '../../contexts/vulcan';
import {getWorkflows} from '../../api/vulcan';
import {defaultInputValues} from '../../utils/workflow';

import SessionManager from './session/session_manager';
import StepsList from './steps/steps_list';

// Hardcode for now, since only one workflow
const WORKFLOW_SHORT = 'umap';
const WORKFLOW_NAME = `${WORKFLOW_SHORT}.cwl`;

export default function WorkflowManager() {
  const invoke = useActionInvoker();
  const {setWorkflow, setPathIndex, setInputs} = useContext(VulcanContext);

  useEffect(() => {
    getWorkflows()
      .then((response) => {
        let currentWorkflow = response.workflows.find(
          (w) => WORKFLOW_NAME === w.name
        );

        setWorkflow(currentWorkflow);

        // Set the default input values
        setInputs(defaultInputValues(currentWorkflow));

        // first path is always the "work" path
        setPathIndex(0);
      })
      .catch((e) => {
        console.error(e);
        invoke(showMessages([e]));
      });
  }, []);

  return (
    <div className='workflow-manager'>
      <div className='workflow-header'>{WORKFLOW_SHORT}</div>
      <div className='step-wrapper'>
        <div className='step-main-pane-wrapper'>
          <SessionManager></SessionManager>
        </div>
        <div className='step-nav-wrapper'>
          <StepsList></StepsList>
        </div>
      </div>
    </div>
  );
}
