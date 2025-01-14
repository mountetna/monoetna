import React, {useState, useCallback} from 'react';

import 'regenerator-runtime/runtime';
import {json_post} from 'etna-js/utils/fetch';

import TableRow from '@material-ui/core/TableRow';
import TableCell from '@material-ui/core/TableCell';
import Typography from '@material-ui/core/Typography';
import Card from '@material-ui/core/Card';
import CardContent from '@material-ui/core/CardContent';
import {makeStyles} from '@material-ui/core/styles';
import CardActions from '@material-ui/core/CardActions';
import SettingsIcon from '@material-ui/icons/Settings';
import LibraryBooksIcon from '@material-ui/icons/LibraryBooksRounded';
import PlayArrowIcon from '@material-ui/icons/PlayArrowRounded';
import DeleteIcon from '@material-ui/icons/Delete';
import LockIcon from '@material-ui/icons/Lock';
import CheckIcon from '@material-ui/icons/Check';
import ErrorIcon from '@material-ui/icons/ErrorOutline';
import ScheduleIcon from '@material-ui/icons/Schedule';
import Grid from '@material-ui/core/Grid';
import Chip from '@material-ui/core/Chip';

import WorkflowButton from './workflow-button';
import RunPane from './run-pane';
import ConfigurePane from './configure-pane';
import RemovePane from './remove-pane';
import LogsPane from './logs-pane';
import SecretsPane from './secrets-pane';
import {formatTime, runTime, runDesc} from './run-state';
import useAsyncWork from 'etna-js/hooks/useAsyncWork';

import {Workflow, Job} from '../polyphemus';

const StatusIcon = ({status}: {status: string}) => {
  let IconComponent: any;
  if (status == 'completed') IconComponent = CheckIcon;
  else if (status == 'message') IconComponent = ErrorIcon;
  else if (status == 'pending') IconComponent = ScheduleIcon;
  else return null;

  return <IconComponent size='small' />;
};

const useStyles = makeStyles((theme) => ({
  title: {
    color: 'goldenrod',
    justifyContent: 'flex-end',
    alignItems: 'center',
    paddingRight: '15px',
    flexBasis: '200px'
  },
  values: {
    flex: '1 1 auto',
    alignItems: 'center',
    width: 'auto'
  },
  heading: {
    borderBottom: '1px solid #ccc'
  },
  workflow: {
    border: '1px solid #ccc',
    marginBottom: '15px',
    height: 'calc(100vh - 160px)'
  },
  workflowrow: {
    marginBottom: '15px',
    cursor: 'pointer',
    '&:hover': {
      background: '#eee'
    }
  },
  savebar: {
    justifyContent: 'space-between'
  },
  buttons: {
    justifyContent: 'space-evenly'
  },
  statusline: {
    borderBottom: '1px solid #eee'
  },
  completed: {
    color: 'green',
    paddingLeft: '0.5rem'
  },
  error: {
    paddingLeft: '0.5rem',
    color: 'red'
  },
  pending: {
    color: 'goldenrod'
  },
  none: {
    color: 'salmon'
  }
}));

export const WorkflowConfigRow = ({
  workflow,
  runtimeConfig,
  status,
  params,
  output,
  run_interval,
  ran_at,
  job,
  onClick
}:{ workflow: Workflow; runtimeConfig: RuntimeConfig; job: Job | undefined; onClick: Function}) => {
  const [mode, setMode] = useState<string | null>(null);
  const [message, setMessage] = useState('');
  const [error, setError] = useState('');

  const { workflow_name } = workflow;

  const classes: any = useStyles();

  return (
    <TableRow className={classes.workflowrow} onClick={ onClick as React.MouseEventHandler<HTMLTableRowElement>}>
      <TableCell>{job?.name}</TableCell>
      <TableCell>{workflow_name}</TableCell>
      <TableCell align="right">{status || 'none'}</TableCell>
      <TableCell align="right">{ran_at ? formatTime(ran_at) : 'never'}</TableCell>
      <TableCell align="right">{runTime(ran_at, run_interval)}</TableCell>
      <TableCell align="right">{runDesc(run_interval)}</TableCell>
    </TableRow>
  );
};

export const WorkflowConfig = ({
  workflow,
  job,
  onUpdate
}:{
  workflow: Workflow;
  job: Job | undefined;
  onUpdate: Function
}) => {
  const [mode, setMode] = useState<string | null>(null);
  const [message, setMessage] = useState('');
  const [error, setError] = useState('');

  const clearMessages = useCallback(() => {
    setMessage('');
    setError('');
  }, []);

  const toggleMode = (m: string) => {
    mode == m ? setMode(null) : setMode(m);
    clearMessages();
  };

  const classes: any = useStyles();

  const [_, postUpdate] = useAsyncWork(
    function postUpdate(update: any) {
      clearMessages();
      return json_post(`/api/workflow/${project_name}/update/${config_id}`, update)
        .then((workflow) => {
          onUpdate(workflow);
          setMessage('Saved!');
        })
        .catch((r) => r.then(({error}: {error: string}) => setError(error)));
    },
    {cancelWhenChange: []}
  );

  return (
    <Card className={classes.workflow} elevation={0} key={workflow}>
      <CardContent>
        <CardActions>
          <Grid direction='row' container>
            <Grid direction='column' container item xs={9}>
              <Grid direction='row' className={classes.statusline} container>
                <Grid container className={classes.title} item>
                  <Typography>Job Type</Typography>
                </Grid>
                <Grid item className={classes.values}>
                  <Typography>
                    {job?.name}
                  </Typography>
                </Grid>
              </Grid>
              <Grid direction='row' className={classes.statusline} container>
                <Grid container className={classes.title} item>
                  <Typography>Last Ran</Typography>
                </Grid>
                <Grid item className={classes.values}>
                  <Typography>
                    {ran_at ? formatTime(ran_at) : 'never'}
                  </Typography>
                </Grid>
              </Grid>
              <Grid direction='row' className={classes.statusline} container>
                <Grid container className={classes.title}>
                  <Typography>Last Status</Typography>
                </Grid>
                <Grid
                  className={`${classes.values} ${classes[status]}`}
                  direction='row'
                  item
                  container
                >
                  <Grid item>
                    <StatusIcon status={status} />
                  </Grid>
                  <Grid item>
                    <Typography>{status || 'none'}</Typography>
                  </Grid>
                </Grid>
              </Grid>
              <Grid direction='row' className={classes.statusline} container>
                <Grid container className={classes.title} item>
                  <Typography>Next Run</Typography>
                </Grid>
                <Grid item className={classes.values}>
                  <Typography>{runTime(ran_at, run_interval)} </Typography>
                </Grid>
              </Grid>
            </Grid>
            <Grid
              item
              container
              direction='row'
              className={classes.buttons}
              alignItems='center'
              xs={3}
            >
              <WorkflowButton selected={mode} mode='run' onClick={toggleMode}>
                <PlayArrowIcon />
              </WorkflowButton>
              <WorkflowButton selected={mode} mode='logs' onClick={toggleMode}>
                <LibraryBooksIcon />
              </WorkflowButton>
              <WorkflowButton selected={mode} mode='configure' onClick={toggleMode}>
                <SettingsIcon />
              </WorkflowButton>
              <WorkflowButton selected={mode} mode='secrets' onClick={toggleMode}>
                <LockIcon />
              </WorkflowButton>
              <WorkflowButton selected={mode} mode='remove' onClick={toggleMode}>
                <DeleteIcon />
              </WorkflowButton>
            </Grid>
          </Grid>
        </CardActions>
        {(error || message) && (
          <Grid item className={error ? classes.error : classes.completed}>
            <Typography>{error || message}</Typography>
          </Grid>
        )}
        <ConfigurePane
          config_id={config_id}
          project_name={project_name}
          selected={mode}
          config={config}
          job={job}
          update={postUpdate}
        />
        <RunPane
          selected={mode}
          run_interval={run_interval}
          update={postUpdate}
          config={config}
          params={params}
          param_opts={job ? job.params : null}
        />
        <RemovePane selected={mode} update={postUpdate} />
        <LogsPane selected={mode} config_id={config_id} name={name} project_name={project_name} />
        <SecretsPane
          selected={mode}
          update={postUpdate}
          keys={job ? job.secrets : null}
          secrets={secrets}
        />
      </CardContent>
    </Card>
  );
};
