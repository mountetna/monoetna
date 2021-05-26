import React, {useCallback, useContext, useEffect, useMemo} from 'react';
import ReactModal from 'react-modal';
import FlatButton from 'etna-js/components/flat-button';

import {VulcanContext} from '../../../contexts/vulcan_context';
import {setSession} from "../../../actions/vulcan";
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
import { readTextFile, downloadBlob } from 'etna-js/utils/blob';

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
  const {state, dispatch, requestPoll} = useContext(VulcanContext);
  const workflow = useContext(WorkflowContext);

  const [modalIsOpen, setIsOpen] = React.useState(false);
  const {steps} = workflow;
  const {status, session} = state;

  const openModal = useCallback(() => setIsOpen(true), [setIsOpen]);
  const closeModal = useCallback(() => setIsOpen(false), [setIsOpen]);

  const saveSession = useCallback(
    () => {
      downloadBlob({
      data: JSON.stringify(session, null, 2),
      filename: `${name}.json`,
      contentType: 'text/json'})
    }, [session]
  );

  const openSession = () => {
    readTextFile('*.json').then(
      session_json => {
        dispatch(setSession(JSON.parse(session_json)))
      }
    )
  }

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
        <FlatButton className='header-btn save' icon='save' label='Save' title='Save workflow parameters to file' onClick={saveSession} disabled={running}/>
        <FlatButton className='header-btn open' icon='folder-open' label='Open' title='Load workflow parameters from file' onClick={openSession} disabled={running}/>
      </div>
      <div className='session-feed-container'>
        <InputFeed />
        <OutputFeed />
      </div>
    </div>
  );
}
