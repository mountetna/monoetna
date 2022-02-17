import React, {
  useCallback,
  useEffect,
  useState,
  useContext,
  useMemo
} from 'react';
import * as _ from 'lodash';
import ReactModal from 'react-modal';
import FlatButton from 'etna-js/components/flat-button';

import {makeStyles} from '@material-ui/core/styles';
import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
import {pushLocation} from 'etna-js/actions/location_actions';
import {showMessages} from 'etna-js/actions/message_actions';

import Breadcrumbs from '@material-ui/core/Breadcrumbs';
import Typography from '@material-ui/core/Typography';
import Grid from '@material-ui/core/Grid';
import TextField from '@material-ui/core/TextField';
import Link from '@material-ui/core/Link';

import {VulcanContext} from '../../../contexts/vulcan_context';
import {
  clearCommittedStepPending,
  setSession,
  setFigure,
  setSessionAndFigure
} from '../../../actions/vulcan_actions';
import InputFeed from './input_feed';
import OutputFeed from './output_feed';
import Vignette from '../vignette';
import {workflowName} from '../../../selectors/workflow_selectors';
import {useWorkflow} from '../../../contexts/workflow_context';
import {readTextFile, downloadBlob} from 'etna-js/utils/blob';
import {defaultSession} from '../../../reducers/vulcan_reducer';
import {VulcanSession} from '../../../api_types';
import {json_post} from 'etna-js/utils/fetch';
import {Debouncer} from 'etna-js/utils/debouncer';
import Tooltip from '@material-ui/core/Tooltip';

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

const useStyles = makeStyles((theme) => ({
  title: {
    width: '800px'
  },
  titleText: {
    textOverflow: 'ellipsis',
    overflow: 'hidden'
  }
}));

export default function SessionManager() {
  const [localTitle, setLocalTitle] = useState('');
  const [saving, setSaving] = useState(false);
  const {state, dispatch, showErrors, requestPoll, cancelPolling} = useContext(
    VulcanContext
  );
  const {workflow, hasPendingEdits, complete} = useWorkflow();

  const [modalIsOpen, setIsOpen] = React.useState(false);
  const {session, figure, committedStepPending} = state;

  const invoke = useActionInvoker();

  const classes = useStyles();

  const name = workflowName(workflow);
  const openModal = useCallback(() => setIsOpen(true), [setIsOpen]);
  const closeModal = useCallback(() => setIsOpen(false), [setIsOpen]);

  const run = useCallback(() => {
    requestPoll(true);
    dispatch(clearCommittedStepPending());
  }, [requestPoll, dispatch]);
  const stop = useCallback(() => cancelPolling(), [cancelPolling]);

  const handleErrorResponse = useCallback(
    (err: Promise<any>) => {
      err.then((e: any) => invoke(showMessages(e.errors || [e.error])));
    },
    [invoke]
  );
  const cancelSaving = () => {
    setSaving(false);
  };

  const saveSession = useCallback(() => {
    if (hasPendingEdits) {
      if (!confirm('Pending edits will be discarded when saving. Proceed?')) {
        return;
      }
    }

    let params = {
      ...figure,
      workflow_name: name,
      inputs: {...session.inputs}
    };

    if (!params.title) {
      params.title = prompt('Set a title for this figure');
      if (!params.title) return;
    }

    setSaving(true);
    if (params.figure_id) {
      json_post(
        `/api/${session.project_name}/figure/${params.figure_id}/update`,
        params
      )
        .then((figureResponse) => {
          setSessionAndFigure(figureResponse);
        })
        .catch(handleErrorResponse)
        .finally(cancelSaving);
    } else {
      json_post(`/api/${session.project_name}/figure/create`, params)
        .then((figure) =>
          invoke(
            pushLocation(`/${figure.project_name}/figure/${figure.figure_id}`)
          )
        )
        .catch(handleErrorResponse)
        .finally(cancelSaving);
    }
  }, [hasPendingEdits, session, name, figure, invoke, handleErrorResponse]);

  const saveSessionToBlob = useCallback(() => {
    if (hasPendingEdits) {
      if (!confirm('Pending edits will be discarded when saving. Proceed?')) {
        return;
      }
    }

    downloadBlob({
      data: JSON.stringify(session, null, 2),
      filename: `${name}.json`,
      contentType: 'text/json'
    });
  }, [hasPendingEdits, session, name]);

  const openSession = () => {
    showErrors(
      readTextFile('*.json').then((sessionJson) => {
        const session: VulcanSession = JSON.parse(sessionJson);
        if (
          session.workflow_name !== state.session.workflow_name ||
          session.project_name !== state.session.project_name
        ) {
          // TODO: Confirm the user has project and workflow access before doing the navigation?
          if (
            confirm(
              'This session file belongs to a different project / workflow combination, navigate to that page?'
            )
          ) {
            location.href =
              location.origin +
              ROUTES.workflow(session.project_name, session.workflow_name);
          }
          throw new Error(
            'Cannot read session file, incompatible with current page.'
          );
        }
        dispatch(setSession(session));
      })
    );
  };

  const resetSession = useCallback(() => {
    const newSession = {
      ...defaultSession,
      workflow_name: session.workflow_name,
      project_name: session.project_name
    };
    dispatch(setSession(newSession));
    requestPoll();
  }, [session, dispatch, requestPoll]);

  const running = state.pollingState > 0;
  const disableRunButton =
    complete || running || (hasPendingEdits && !committedStepPending);

  useEffect(() => {
    if (figure.title) setLocalTitle(figure.title);
  }, [figure]);

  const [debouncer, setDebouncer] = useState(new Debouncer({windowMs: 800}));

  const debouncedSetTitle = useCallback(
    (newTitle: any) => {
      setLocalTitle(newTitle);
      debouncer.ready(() => {
        dispatch(setFigure({...figure, title: newTitle}));
      });
    },
    [figure, debouncer, dispatch]
  );

  const inputsChanged = useMemo(() => {
    return !_.isEqual(figure.inputs, session.inputs);
  }, [figure, session]);

  const canSave = useMemo(() => {
    return inputsChanged && !(running || saving);
  }, [running, saving, inputsChanged]);

  if (!name || !session) return null;

  return (
    <div className='session-manager'>
      <div className='session-header'>
        <Breadcrumbs className='session-workflow-name'>
          <Link href={`/${session.project_name}`}>{session.project_name}</Link>
          <Typography>{workflow.displayName}</Typography>
          <Grid container className={classes.title}>
            <Tooltip title={localTitle}>
              <TextField
                fullWidth
                value={localTitle}
                margin='none'
                InputProps={{
                  disableUnderline: true,
                  inputProps: {
                    className: classes.titleText
                  }
                }}
                variant='standard'
                onChange={(e) => debouncedSetTitle(e.target.value)}
                placeholder='Untitled'
              />
            </Tooltip>
          </Grid>
        </Breadcrumbs>
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
        {state.pollingState ? (
          <FlatButton
            className={`header-btn`}
            icon='stop'
            label='Stop'
            title='Stop workflow'
            onClick={stop}
          />
        ) : (
          <FlatButton
            className={`header-btn run`}
            icon='play'
            label='Run'
            title='Run workflow'
            onClick={run}
            disabled={disableRunButton}
          />
        )}
        <FlatButton
          className='header-btn save'
          icon='save'
          label='Save'
          title='Save current workflow parameters to current figure'
          onClick={saveSession}
          disabled={!canSave}
        />
      </div>
      <div className='session-feed-container'>
        <InputFeed />
        <OutputFeed />
      </div>
    </div>
  );
}
