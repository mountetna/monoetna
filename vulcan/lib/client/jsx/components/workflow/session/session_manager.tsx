import React, {useCallback, useContext, useEffect, useMemo} from 'react';
import ReactModal from 'react-modal';
import FlatButton from 'etna-js/components/flat-button';

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
  const context = useContext(VulcanContext);
  const {state, requestPoll} = context;

  const workflow = state.workflow;
  if (!workflow) return null;
  const name = workflowName(workflow);
  if (!name) return null;

  const [modalIsOpen, setIsOpen] = React.useState(false);
  function openModal() {
    setIsOpen(true);
  }

  function closeModal() {
    setIsOpen(false);
  }

  const {steps} = workflow;
  const {status, session} = state;

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

  return (
    <div className='session-manager'>
      <div className='session-header'>
        <span className='session-workflow-name'>
          {workflow.description || name} - {session.project_name}
        </span>
        {workflow.vignette && (
          <React.Fragment>
            <FlatButton icon='book' className='header-btn vignette' label='Vignette' onClick={openModal}/>
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
        <FlatButton className='header-btn run' icon='play' label='Run' title='Run workflow' onClick={run} disabled={disableRunButton} />
        <FlatButton className='header-btn save' icon='save' label='Save' title='Save workflow parameters to file' onClick={saveSession}/>
        <FlatButton className='header-btn open' icon='folder-open' label='Open' title='Load workflow parameters from file' onClick={openSession}/>
      </div>
      <div className='session-feed-container'>
        <InputFeed />
        <OutputFeed />
      </div>
    </div>
  );
}
