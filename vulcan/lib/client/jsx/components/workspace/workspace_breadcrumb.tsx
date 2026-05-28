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

import Breadcrumbs from '@material-ui/core/Breadcrumbs';
import Typography from '@material-ui/core/Typography';
import TextField from '@material-ui/core/TextField';
import Link from '@material-ui/core/Link';
import Tooltip from '@material-ui/core/Tooltip';

import {WorkspaceContext} from '../../contexts/workspace_context';
import {
  clearRunTriggers,
  setAttemptingToRun,
  setRunning,
  setWorkspace,
} from '../../actions/vulcan_actions';
import InputFeed from './input_feed';
import OutputFeed from './output_feed';
// import Vignette from '../vignette';
import VulcanHelp from './drawers/vulcan_help';
import { hasScheduledSteps, workflowName, workspaceFromRaw } from '../../selectors/workflow_selectors';
import {useWorkspace} from '../../contexts/workspace_context';
// import {json_get} from 'etna-js/utils/fetch';
import useUserHooks from '../../contexts/useUserHooks';
import Tag from '../dashboard/tag';
import Grid from '@material-ui/core/Grid';
import { useDataSync, useRunSyncing } from './data_sync';
import Vignette from './drawers/vignette';

// import RevisionHistory from 'etna-js/components/revision-history';

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
  },
  tags: {
    padding: '12.5px !important'
  },
  editTags: {
    width: '600px'
  }
}));

const WorkspaceBreadcrumb = ({editor}) => {
  const classes = useStyles();

  const { state: { workflow, workspace } } = React.useContext(WorkspaceContext);

  return (
    <Breadcrumbs
      className='session-workflow-name'
      classes={{
        li: classes.title
      }}
    >
      <Link href={`/${workflow.project_name}`}>{workflow.project_name}</Link>
      <Typography>{workflow.name}</Typography>
      {editor ? (
        <Grid container direction='row'>
          <Grid item>
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
              disabled={updating}
              variant='standard'
              onChange={(e) => setLocalTitle(e.target.value)}
              placeholder='Untitled'
            />
          </Grid>
          {localTitle != workspace.name && <Grid item xs={2}>
            <FlatButton
              className='header-btn-name-save'
              icon={updatingTitle? 'spinner fa-spin' : 'save'}
              label='Save'
              title='Save Workspace Title'
              disabled={updating}
              onClick={() => {
                setUpdatingTitle(true);
                handleUpdateWorkspace(localTitle, undefined)
              }}
            />
          </Grid>}
        </Grid>
      ) : (
        <Typography>{workspace.name}</Typography>
      )
      }
    </Breadcrumbs>
  );
}

export default WorkspaceBreadcrumb;
