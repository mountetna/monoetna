import React, {useState, useCallback} from 'react';
import {json_post} from 'etna-js/utils/fetch';

import Typography from '@material-ui/core/Typography';
import Card from '@material-ui/core/Card';
import CardContent from '@material-ui/core/CardContent';
import {makeStyles} from '@material-ui/core/styles';
import CardActions from '@material-ui/core/CardActions';
import ButtonGroup from '@material-ui/core/ButtonGroup';
import Button from '@material-ui/core/Button';
import SettingsIcon from '@material-ui/icons/Settings';
import LibraryBooksIcon from '@material-ui/icons/LibraryBooksRounded';
import PlayArrowIcon from '@material-ui/icons/PlayArrowRounded';
import DeleteIcon from '@material-ui/icons/Delete';
import LockIcon from '@material-ui/icons/Lock';
import ExpandMoreIcon from '@material-ui/icons/ExpandMoreRounded';
import CheckIcon from '@material-ui/icons/Check';
import ErrorIcon from '@material-ui/icons/ErrorOutline';
import ScheduleIcon from '@material-ui/icons/Schedule';
import Box from '@material-ui/core/Box';
import Select from '@material-ui/core/Select';
import TextField from '@material-ui/core/TextField';
import InputAdornment from '@material-ui/core/InputAdornment';
import MenuItem from '@material-ui/core/MenuItem';
import Grid from '@material-ui/core/Grid';
import Collapse from '@material-ui/core/Collapse';

import EtlButton from './etl-button';
import RunPane from './run-pane';
import ConfigurePane from './configure-pane';
import RemovePane from './remove-pane';
import LogsPane from './logs-pane';
import SecretsPane from './secrets-pane';
import {formatTime, runTime} from './run-state';

const StatusIcon = ({status}: {status: string}) => {
  let IconComponent: any;
  if (status == 'completed') IconComponent = CheckIcon;
  else if (status == 'error') IconComponent = ErrorIcon;
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
  etl: {
    border: '1px solid #ccc',
    marginBottom: '15px'
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
    color: 'green'
  },
  error: {
    color: 'red'
  },
  pending: {
    color: 'goldenrod'
  },
  none: {
    color: 'salmon'
  }
}));

const EtlConfig = ({
  project_name,
  etl,
  name,
  config,
  status,
  secrets,
  params,
  output,
  run_interval,
  ran_at,
  job,
  onUpdate
}: Etl & {job: Job | undefined; onUpdate: Function}) => {
  const [mode, setMode] = useState<string | null>(null);
  const [error, setError] = useState('');
  const toggleMode = (m: string) => (mode == m ? setMode(null) : setMode(m));

  const classes: any = useStyles();

  const postUpdate = useCallback(
    (update: any) => {
      return json_post(`/api/etl/${project_name}/update/${name}`, update)
        .then((etl) => {
          setMode(null);
          onUpdate(etl);
          setError('');
        })
        .catch((r) => r.then(({error}: {error: string}) => setError(error)));
    },
    [project_name, name]
  );

  return (
    <Card className={classes.etl} elevation={0} key={etl}>
      <CardContent>
        <Typography className={classes.heading}>{name}</Typography>

        <CardActions>
          <Grid direction='row' container>
            <Grid direction='column' container item xs={9}>
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
              <EtlButton selected={mode} mode='run' onClick={toggleMode}>
                <PlayArrowIcon />
              </EtlButton>
              <EtlButton selected={mode} mode='logs' onClick={toggleMode}>
                <LibraryBooksIcon />
              </EtlButton>
              <EtlButton selected={mode} mode='configure' onClick={toggleMode}>
                <SettingsIcon />
              </EtlButton>
              <EtlButton selected={mode} mode='secrets' onClick={toggleMode}>
                <LockIcon />
              </EtlButton>
              <EtlButton selected={mode} mode='remove' onClick={toggleMode}>
                <DeleteIcon />
              </EtlButton>
            </Grid>
          </Grid>
        </CardActions>
        {error && (
          <Grid item className={classes.error}>
            <Typography>{error}</Typography>
          </Grid>
        )}
        <ConfigurePane
          name={name}
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
        <LogsPane selected={mode} name={name} project_name={project_name} />
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

export default EtlConfig;
