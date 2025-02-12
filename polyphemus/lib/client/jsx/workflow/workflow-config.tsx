import React, {useState, useCallback} from 'react';

import 'regenerator-runtime/runtime';
import {json_post} from 'etna-js/utils/fetch';

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

import {Workflow, Status, Job, Runtime} from '../polyphemus';

const StatusIcon = ({status}: {status: string}) => {
  let IconComponent: any;
  if (status == 'succeeded') IconComponent = CheckIcon;
  else if (status == 'failed') IconComponent = ErrorIcon;
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

export const WorkflowConfig = ({
  workflow,
  status,
  runtime,
  job,
  onUpdateWorkflow,
  onUpdateRuntime,
}:{
  workflow: Workflow;
  runtime: Runtime;
  status: Status;
  job: Job | undefined;
  onUpdateWorkflow: Function
  onUpdateRuntime: Function
}) => {
  const [mode, setMode] = useState<string | null>(null);
  const [message, setMessage] = useState('');
  const [error, setError] = useState('');

  const { project_name, config_id, config, secrets } = workflow;

  const { run_interval, config: params={} } = runtime;

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
      return json_post(`/api/workflows/${project_name}/update/${config_id}`, update)
        .then((workflow) => {
          onUpdateWorkflow(workflow);
          setMessage('Saved!');
        })
        .catch((r) => r.then(({error}: {error: string}) => setError(error)));
    },
    {cancelWhenChange: []}
  );

  const [__, postRuntimeUpdate] = useAsyncWork(
    function postRuntimeUpdate(update: any) {
      clearMessages();
      return json_post(`/api/workflows/${project_name}/runtime_configs/update/${config_id}`, update)
        .then((runtime) => {
          onUpdateRuntime(runtime);
          setMessage('Saved!');
        })
        .catch((r) => r.then(({error}: {error: string}) => setError(error)));
    },
    {cancelWhenChange: []}
  );

  return (
    <Card className={classes.workflow} elevation={0} key={workflow.config_id}>
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
                  <Typography>Last Completed</Typography>
                </Grid>
                <Grid item className={classes.values}>
                  <Typography>
                    { status.pipeline_finished_at || 'never run'}
                  </Typography>
                </Grid>
              </Grid>
              <Grid direction='row' className={classes.statusline} container>
                <Grid container className={classes.title}>
                  <Typography>Last Status</Typography>
                </Grid>
                <Grid
                  className={`${classes.values} ${classes[status.pipeline_state]}`}
                  direction='row'
                  item
                  container
                >
                  <Grid item>
                    <StatusIcon status={status.pipeline_state} />
                  </Grid>
                  <Grid item>
                    <Typography>{status.pipeline_state || 'never run'}</Typography>
                  </Grid>
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
          config_id={config_id}
          project_name={project_name}
          run_interval={run_interval}
          update={postRuntimeUpdate}
          config={config}
          params={params}
          param_opts={job ? job.runtime_params : null}
        />
        <RemovePane selected={mode} update={postRuntimeUpdate} />
        <LogsPane selected={mode} config_id={config_id} project_name={project_name} />
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
