import React, {useCallback, useContext, useEffect, useMemo} from 'react';
import ReactModal from 'react-modal';
import Icon from 'etna-js/components/icon';

import {VulcanContext} from '../../../contexts/vulcan_context';
import InputFeed from './input_feed';
import OutputFeed from './output_feed';
import Vignette from '../vignette';
import {
  allWorkflowPrimaryInputSources,
  statusOfStep,
  uiOutputOfStep,
  workflowName
} from '../../../selectors/workflow_selectors';
import {WorkflowContext} from "../../../contexts/workflow_context";

const modalStyles = {
  content: {
    top: '50%',
    left: '50%',
    right: 'auto',
    bottom: 'auto',
    marginRight: '-50%',
    transform: 'translate(-50%, -50%)'
  }
};

export default function SessionManager() {
  const {state, requestPoll} = useContext(VulcanContext);
  const workflow = useContext(WorkflowContext);
  const [modalIsOpen, setIsOpen] = React.useState(false);
  const {steps} = workflow;
  const {status, session} = state;

  const openModal = useCallback(() => setIsOpen(true), [setIsOpen]);
  const closeModal = useCallback(() => setIsOpen(false), [setIsOpen]);

  // We are done once every step either has a download or that step is a uiOutput.
  const complete = useMemo(
    () =>
      steps[0].every(
        (step) => uiOutputOfStep(step) || statusOfStep(step, status)?.downloads
      ),
    [steps, status]
  );

  const idle = useMemo(
    () =>
      steps[0].every(
        (step) =>
          uiOutputOfStep(step) ||
          statusOfStep(step, status)?.status !== 'running'
      ),
    [steps, status]
  );

  const primaryInputsReady = useMemo(
    () =>
      allWorkflowPrimaryInputSources(workflow).every(
        (source) => source in state.inputs
      ),
    [state.inputs, workflow]
  );

  const running = !idle;

  const run = useCallback(() => {
    requestPoll(true);
  }, [requestPoll]);

  const disableRunButton = complete || running || !primaryInputsReady;

  const name = workflowName(workflow);
  if (!name) return null;

  return (
    <div className='session-manager'>
      <div className='session-header'>
        <span className='session-workflow-name'>
          {workflow.description || name} - {session.project_name}
        </span>
        {workflow.vignette && (
          <React.Fragment>
            <div className='header-btn' onClick={openModal}>
              <div className='vignette-btn'>
                Vignette
                <Icon className='vignette' icon='book' />
              </div>
            </div>
            <ReactModal
              isOpen={modalIsOpen}
              onRequestClose={closeModal}
              style={modalStyles}
              contentLabel='Vignette'
            >
              <Vignette workflowName={name} />
            </ReactModal>
          </React.Fragment>
        )}
        <div
          onClick={run}
          className={`run-workflow-btn ${
            disableRunButton ? 'disabled' : ''
          } header-btn`}
        >
          Run
          <Icon
            className='run'
            disabled={complete || running || !primaryInputsReady}
            title='Run workflow'
            icon='play'
          />
        </div>
      </div>
      <div className='session-feed-container'>
        <InputFeed />
        <OutputFeed />
      </div>
    </div>
  );
}
