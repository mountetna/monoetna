import React, {useCallback, useContext} from 'react';
import ReactModal from 'react-modal';
import FlatButton from 'etna-js/components/flat-button';

import {VulcanContext} from '../../../contexts/vulcan_context';
import {clearCommittedStepPending, setSession} from '../../../actions/vulcan_actions';
import InputFeed from './input_feed';
import OutputFeed from './output_feed';
import Vignette from '../vignette';
import {
  workflowName
} from '../../../selectors/workflow_selectors';
import {useWorkflow} from '../../../contexts/workflow_context';
import {readTextFile, downloadBlob} from 'etna-js/utils/blob';
import { defaultSession } from '../../../reducers/vulcan_reducer';
import {VulcanSession} from "../../../api_types";

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
  const {state, dispatch, showErrors, requestPoll, cancelPolling} = useContext(VulcanContext);
  const {workflow, hasPendingEdits, complete} = useWorkflow();

  const [modalIsOpen, setIsOpen] = React.useState(false);
  const {session, committedStepPending} = state;

  const name = workflowName(workflow);
  const openModal = useCallback(() => setIsOpen(true), [setIsOpen]);
  const closeModal = useCallback(() => setIsOpen(false), [setIsOpen]);

  const run = useCallback(() => {requestPoll(true); dispatch(clearCommittedStepPending())}, [requestPoll]);
  const stop = useCallback(() => cancelPolling(), [cancelPolling]);

  const saveSession = useCallback(() => {
    if (hasPendingEdits) {
      alert('You have uncommitted changes, reset or commit those before saving.')
      return;
    }

    downloadBlob({
      data: JSON.stringify(session, null, 2),
      filename: `${name}.json`,
      contentType: 'text/json'
    });
  }, [hasPendingEdits, session, name]);

  const openSession = () => {
    showErrors(readTextFile('*.json').then((sessionJson) => {
      const session: VulcanSession = JSON.parse(sessionJson);
      if (session.workflow_name !== state.session.project_name || session.project_name !== state.session.project_name) {
        // TODO: Confirm the user has project and workflow access before doing the navigation?
        if (confirm('This session file belongs to a different project / workflow combination, navigate to that page?')) {
          location.href = location.origin + ROUTES.workflow(session.project_name, session.workflow_name);
        }
        throw new Error('Cannot read session file, incompatible with current page.')
      }
      dispatch(setSession(session));
    }));
  };

  const resetSession = useCallback(
    () => {
      const newSession = {
        ...defaultSession,
        workflow_name: session.workflow_name,
        project_name: session.project_name,
      }
      dispatch(setSession(newSession));
      requestPoll();
    }, [session.workflow_name, session.project_name, dispatch, requestPoll]
  );

  const running = state.pollingState > 0;
  const disableRunButton = complete || running || (hasPendingEdits && !committedStepPending);

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
        { state.pollingState ? <FlatButton
          className={`header-btn`}
          icon='stop'
          label='Stop'
          title='Stop workflow'
          onClick={stop}
        /> :
        <FlatButton
          className={`header-btn run`}
          icon='play'
          label='Run'
          title='Run workflow'
          onClick={run}
          disabled={disableRunButton}
        /> }
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
        <FlatButton
          className='header-btn reset'
          icon='undo'
          label='Reset'
          title='Reset to the default workflow parameters'
          onClick={resetSession}
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
