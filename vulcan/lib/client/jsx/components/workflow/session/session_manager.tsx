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

import Breadcrumbs from '@material-ui/core/Breadcrumbs';
import Typography from '@material-ui/core/Typography';
import TextField from '@material-ui/core/TextField';
import Link from '@material-ui/core/Link';
import Tooltip from '@material-ui/core/Tooltip';

import {VulcanContext} from '../../../contexts/vulcan_context';
import {
  clearCommittedStepPending,
  setSession,
  setSessionAndFigure
} from '../../../actions/vulcan_actions';
import InputFeed from './input_feed';
import OutputFeed from './output_feed';
import Vignette from '../vignette';
import {workflowName} from '../../../selectors/workflow_selectors';
import {useWorkflow} from '../../../contexts/workflow_context';
import {readTextFile, downloadBlob} from 'etna-js/utils/blob';
import {defaultSession} from '../../../reducers/vulcan_reducer';
import {
  VulcanFigure,
  VulcanFigureSession,
  VulcanSession
} from '../../../api_types';
import useUserHooks from '../../useUserHooks';

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
    '&:last-child': {flex: '1 1 auto'}
  },
  titleText: {
    textOverflow: 'ellipsis',
    overflow: 'hidden'
  }
}));

export default function SessionManager() {
  const [saving, setSaving] = useState(false);
  const {
    state,
    dispatch,
    showErrors,
    requestPoll,
    cancelPolling,
    updateFigure,
    createFigure,
    clearLocalSession
  } = useContext(VulcanContext);
  const {workflow, hasPendingEdits, complete} = useWorkflow();
  const {canEdit} = useUserHooks();

  const [modalIsOpen, setIsOpen] = React.useState(false);
  const {session, figure, committedStepPending} = state;

  const [tags, setTags] = useState<string[]>(figure.tags || []);
  const [localTitle, setLocalTitle] = useState(figure.title);

  const invoke = useActionInvoker();

  const classes = useStyles();

  const name = workflowName(workflow);
  const openModal = useCallback(() => setIsOpen(true), [setIsOpen]);
  const closeModal = useCallback(() => setIsOpen(false), [setIsOpen]);

  const run = useCallback(() => {
    showErrors(requestPoll(true));
    dispatch(clearCommittedStepPending());
  }, [requestPoll, dispatch, showErrors]);
  const stop = useCallback(() => cancelPolling(), [cancelPolling]);

  const cancelSaving = () => {
    setSaving(false);
  };

  const handleSaveOrCreate = useCallback(
    (figure: VulcanFigure) => {
      let params = {
        ...figure,
        workflow_name: name,
        inputs: {...session.inputs},
        title: localTitle,
        tags: [...tags]
      };

      if (!params.title) {
        params.title = prompt('Set a title for this figure');
        if (!params.title) return;
      }

      setSaving(true);
      if (params.figure_id) {
        showErrors(
          updateFigure(session.project_name, params)
            .then((figureResponse) => {
              dispatch(setSessionAndFigure(figureResponse));
            })
            .finally(cancelSaving)
        );
      } else {
        showErrors(
          createFigure(session.project_name, params).then(
            (figure: VulcanFigureSession) => {
              cancelSaving();
              clearLocalSession(
                figure.workflow_name,
                figure.project_name,
                null
              );
              invoke(
                pushLocation(
                  `/${figure.project_name}/figure/${figure.figure_id}`
                )
              );
            }
          )
        );
      }
    },
    [
      name,
      localTitle,
      session,
      tags,
      cancelSaving,
      showErrors,
      updateFigure,
      invoke,
      dispatch,
      clearLocalSession,
      pushLocation
    ]
  );

  const saveSession = useCallback(() => {
    if (hasPendingEdits) {
      if (!confirm('Pending edits will be discarded when saving. Proceed?')) {
        return;
      }
    }

    handleSaveOrCreate(figure);
  }, [hasPendingEdits, handleSaveOrCreate, figure]);

  const copyFigure = useCallback(() => {
    let clone = {
      ...figure
    };

    delete clone.figure_id;

    handleSaveOrCreate(clone);
  }, [figure, handleSaveOrCreate]);

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

  const inputsChanged = useMemo(() => {
    return !_.isEqual(figure.inputs, session.inputs);
  }, [figure, session]);

  const titleChanged = useMemo(() => {
    return localTitle !== figure.title;
  }, [figure, localTitle]);

  const tagsChanged = useMemo(() => {
    return !_.isEqual(tags, figure.tags);
  }, [figure, tags]);

  const canSave = useMemo(() => {
    return (
      (titleChanged || inputsChanged || tagsChanged) && !(running || saving)
    );
  }, [running, saving, inputsChanged, titleChanged, tagsChanged]);

  const editor = useMemo(() => canEdit(figure) || !figure.figure_id, [figure]);

  const isPublic = useMemo(() => (tags || []).includes('public'), [tags]);

  if (!name || !session) return null;
  console.log('tags', tags, isPublic, figure);
  return (
    <div className='session-manager'>
      <div className='session-header'>
        <Breadcrumbs
          className='session-workflow-name'
          classes={{
            li: classes.title
          }}
        >
          <Link href={`/${session.project_name}`}>{session.project_name}</Link>
          <Typography>{workflow.displayName}</Typography>
          <Tooltip title={localTitle || ''}>
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
              onChange={(e) => setLocalTitle(e.target.value)}
              placeholder='Untitled'
            />
          </Tooltip>
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
        {editor ? (
          <FlatButton
            className='header-btn save'
            icon='save'
            label='Save'
            title='Save current workflow parameters to current figure'
            onClick={saveSession}
            disabled={!canSave}
          />
        ) : (
          <FlatButton
            className='header-btn copy'
            icon='copy'
            label='Copy'
            title='Copy current workflow parameters to new figure'
            onClick={copyFigure}
            disabled={!canSave}
          />
        )}
        {editor ? (
          <FlatButton
            className='header-btn public-private'
            icon={`${isPublic ? 'lock' : 'unlock'}`}
            label={`Make ${isPublic ? 'private' : 'public'}`}
            title={`Make the current figure ${isPublic ? 'private' : 'public'}`}
            onClick={() => setTags(isPublic ? [] : ['public'])}
          />
        ) : null}
      </div>
      <div className='session-feed-container'>
        <InputFeed />
        <OutputFeed />
      </div>
    </div>
  );
}
