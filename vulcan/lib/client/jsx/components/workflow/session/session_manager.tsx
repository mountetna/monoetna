import React, {useCallback, useContext, useEffect, useMemo} from 'react';
import ReactModal from 'react-modal';
import FlatButton from 'etna-js/components/flat-button';
import {showMessages} from 'etna-js/actions/message_actions';

import {VulcanContext} from '../../../contexts/vulcan_context';
import {setSession} from '../../../actions/vulcan';
import InputFeed from './input_feed';
import OutputFeed from './output_feed';
import Vignette from '../vignette';
import {
  allWorkflowPrimaryInputSources,
  statusOfStep,
  uiOutputOfStep,
  workflowName
} from '../../../selectors/workflow_selectors';
import {useWorkflow} from '../../../contexts/workflow_context';
import {readTextFile, downloadBlob} from 'etna-js/utils/blob';

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
  const {state, dispatch, requestPoll, useActionInvoker} = useContext(
    VulcanContext
  );
  const workflow = useWorkflow();
  const invoke = useActionInvoker();

  const [modalIsOpen, setIsOpen] = React.useState(false);
  const {steps} = workflow;
  const {status, session, inputs} = state;

  const name = workflowName(workflow);

  const openModal = useCallback(() => setIsOpen(true), [setIsOpen]);
  const closeModal = useCallback(() => setIsOpen(false), [setIsOpen]);

  const saveSession = useCallback(() => {
    downloadBlob({
      data: JSON.stringify(session, null, 2),
      filename: `${name}.json`,
      contentType: 'text/json'
    });
  }, [session, name]);

  const openSession = () => {
    readTextFile('*.json').then((sessionJson) => {
      dispatch(setSession(JSON.parse(sessionJson)));
    });
  };

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

  const hasValidationErrors = useMemo(
    () => Object.keys(state.validationErrors).length > 0,
    [state.validationErrors]
  );

  const running = !idle;

  const run = useCallback(() => {
    if (hasValidationErrors) {
      invoke(
        showMessages(
          Object.entries(state.validationErrors)
            .map(([inputName, errors]: [string, string[]]) =>
              errors.map((e: string) => `${inputName}: ${e}`)
            )
            .flat()
        )
      );
    } else {
      requestPoll(true);
    }
  }, [requestPoll, hasValidationErrors, invoke, state.validationErrors]);

  const disableRunButton = complete || running || !primaryInputsReady;

  if (!name) return null;

  return (
    <div className='session-manager'>
      <div className='session-header'>
        <span className='session-workflow-name'>
          {workflow.description || name} - {session.project_name}
        </span>
        {workflow.vignette && (
          <React.Fragment>
            <FlatButton
              icon='book'
              className='header-btn vignette'
              label='Vignette'
              onClick={openModal}
            />
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
        <FlatButton
          className='header-btn run'
          icon='play'
          label='Run'
          title='Run workflow'
          onClick={run}
          disabled={disableRunButton}
        />
        <FlatButton
          className='header-btn save'
          icon='save'
          label='Save'
          title='Save workflow parameters to file'
          onClick={saveSession}
          disabled={running}
        />
        <FlatButton
          className='header-btn open'
          icon='folder-open'
          label='Open'
          title='Load workflow parameters from file'
          onClick={openSession}
          disabled={running}
        />
      </div>
      <div className='session-feed-container'>
        <InputFeed />
        <OutputFeed />
      </div>
    </div>
  );
}
