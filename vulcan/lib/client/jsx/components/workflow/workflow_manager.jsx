import React, {useContext, useEffect} from 'react';
import Modal from 'react-modal';

import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
import {showMessages} from 'etna-js/actions/message_actions';

import {VulcanContext} from '../../contexts/vulcan';
import {getWorkflows} from '../../api/vulcan';
import {flatten} from '../../utils/workflow';

import SessionManager from './session/session_manager';
import StepsList from './steps/steps_list';

// Hardcode for now, since only one workflow
const WORKFLOW_SHORT = 'umap';
const WORKFLOW_NAME = `${WORKFLOW_SHORT}.cwl`;

export default function WorkflowManager() {
  const invoke = useActionInvoker();
  const {calculating, setWorkflow, setPathIndex} = useContext(VulcanContext);

  Modal.setAppElement('#root');

  const customStyles = {
    content: {
      top: '50%',
      left: '50%',
      right: 'auto',
      bottom: 'auto',
      marginRight: '-50%',
      transform: 'translate(-50%, -50%)'
    }
  };

  useEffect(() => {
    getWorkflows()
      .then((response) => {
        let currentWorkflow = response.workflows.find(
          (w) => WORKFLOW_NAME === w.name
        );

        // TODO: REMOVE Flatten the workflow steps
        currentWorkflow.steps.push(flatten(currentWorkflow.steps));

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
      <div className='workflow-header'>{WORKFLOW_SHORT}</div>
      <div className='step-wrapper'>
        <div className='step-main-pane-wrapper'>
          <SessionManager></SessionManager>
        </div>
        <div className='step-nav-wrapper'>
          <StepsList></StepsList>
        </div>
      </div>
      <Modal
        isOpen={calculating}
        style={customStyles}
        contentLabel='Calculating Notice'
      >
        <h2>Calculations in progress</h2>
        <div>Archimedes calculations in progress. Please be patient.</div>
      </Modal>
    </div>
  );
}
