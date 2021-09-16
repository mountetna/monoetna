import React, { useState, useCallback } from 'react';
import { json_post } from 'etna-js/utils/fetch';

import ConfigScript from './config-script';

import Typography from '@material-ui/core/Typography';
import Card from '@material-ui/core/Card';
import CardContent from '@material-ui/core/CardContent';
import {makeStyles} from '@material-ui/core/styles';
import CardActions from '@material-ui/core/CardActions';
import ButtonGroup from '@material-ui/core/ButtonGroup';
import Button from '@material-ui/core/Button';
import IconButton from '@material-ui/core/IconButton';
import SettingsIcon from '@material-ui/icons/Settings';
import HistoryIcon from '@material-ui/icons/HistoryRounded';
import LibraryBooksIcon from '@material-ui/icons/LibraryBooksRounded';
import PlayArrowIcon from '@material-ui/icons/PlayArrowRounded';
import DeleteIcon from '@material-ui/icons/Delete';
import ExpandMoreIcon from '@material-ui/icons/ExpandMoreRounded';
import CheckIcon from '@material-ui/icons/Check';
import ErrorIcon from '@material-ui/icons/Error';
import ScheduleIcon from '@material-ui/icons/Schedule';
import Tooltip from '@material-ui/core/Typography';
import Box from '@material-ui/core/Box';
import Select from '@material-ui/core/Select';
import TextField from '@material-ui/core/TextField';
import InputAdornment from '@material-ui/core/InputAdornment';
import MenuItem from '@material-ui/core/MenuItem';
import Grid from '@material-ui/core/Grid';
import Collapse from '@material-ui/core/Collapse';
import Accordion from '@material-ui/core/Accordion';
import AccordionSummary from '@material-ui/core/AccordionSummary';
import AccordionDetails from '@material-ui/core/AccordionDetails';

import { Controlled } from 'react-codemirror2';

const StatusIcon = ({status}:{status:string}) => {
  let IconComponent:any;
  if (status == 'completed') IconComponent = CheckIcon;
  else if (status == 'error') IconComponent = ErrorIcon;
  else if (status == 'pending') IconComponent = ScheduleIcon;
  else IconComponent = ErrorIcon;
  
  return <IconComponent size='small'/>
}

const RUN_ONCE = 0;
const RUN_NEVER = -1;
const RUN_INTERVAL = 1;

const getRunState = (run) => (run != RUN_ONCE && run != RUN_NEVER ? RUN_INTERVAL : run);
const getRunIntervalTime = (run) => (run != RUN_ONCE && run != RUN_NEVER ? run : 600);
const formatTime = (time) => Intl.DateTimeFormat('en-US', {
   year: 'numeric', month: 'long', day: 'numeric',
   hour: 'numeric', minute: 'numeric'
}).format(new Date(time))
const runTime = (ran_at, run) => {
  if (run == RUN_NEVER) return 'never';
  if (run == RUN_ONCE || !ran_at) return 'pending';

  const ranAt = new Date(ran_at);
  const runAt = new Date( ranAt.getTime() + run * 1000 );
  const now = new Date();

  if (runAt < now) return 'pending';

  return formatTime( runAt );
}

const SetRunState = ({run_interval,update}:{run_interval:number, update:Function}) => {
  const [ runState, setRunState ] = useState(getRunState(run_interval));
  const [ runIntervalTime, setRunIntervalTime ] = useState(getRunIntervalTime(run_interval));

  const runValue = () => runState == RUN_INTERVAL ? runIntervalTime : runState;
  const reset = () => { setRunState(getRunState(run_interval)); setRunIntervalTime(getRunIntervalTime(run_interval)); };

  return <Grid spacing={1} container alignItems='center'>
    <Grid item>Run</Grid>
    <Grid item><Select labelId='label' id='select' value={runState} onChange={e => setRunState(parseInt(e.target.value as string))}>
      <MenuItem value={RUN_ONCE}>Once</MenuItem>
      <MenuItem value={RUN_NEVER}>Never</MenuItem>
      <MenuItem value={RUN_INTERVAL}>Interval</MenuItem>
    </Select>
    </Grid>
    {
      runState == RUN_INTERVAL && <Grid item> <TextField
        size='small'
        InputProps={{
          endAdornment:<InputAdornment position='end'>seconds</InputAdornment>
        }}
        value={runIntervalTime}
        onChange={ e => { setRunIntervalTime(parseInt(e.target.value as string)); } }
        />
      </Grid>
    }
    {
      runState != run_interval && (runState != RUN_INTERVAL || runIntervalTime != run_interval) && <React.Fragment>
        <Grid item><Button onClick={ () => update({run_interval: runValue()}) }>Save</Button></Grid>
        <Grid item><Button onClick={ reset } color='secondary'>Reset</Button></Grid>
      </React.Fragment>
    }
  </Grid>
}

const EtlButton = ({children, mode, onClick}:{children:React.ReactNode, mode:string, onClick:Function}) => (
  <Tooltip title={mode}>
    <IconButton onClick={() => onClick(mode)} size='small' aria-label={mode}>
      {children}
    </IconButton>
  </Tooltip>
);

const useStyles = makeStyles( theme => ({
  title: {
    color: 'goldenrod',
    justifyContent: 'flex-end',
    paddingRight: '15px'
  },
  heading: {
    borderBottom: '1px solid #ccc'
  },
  etl: {
    border: '1px solid #ccc',
    marginBottom: '15px',
  },
  buttons: {
    justifyContent: 'space-evenly'
  },
  editor: {
    border: '1px solid #ccc'
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
  }
}));

const EtlConfig = ({project_name, etl,name,config,status,output,run_interval,ran_at,schema,onUpdate}: Etl & {schema:any,onUpdate:Function}) => {
  const [ mode, setMode ] = useState<string | null>(null);
  const toggleMode = (m:string) => mode == m ? setMode(null) : setMode(m);

  const classes:any = useStyles();

  const postUpdate = update => json_post(
    `/api/etl/${project_name}/update/${name}`, update
  ).then( etl => onUpdate(etl) );

  return <Card className={classes.etl} elevation={0} key={etl}>
    <CardContent>
      <Typography className={classes.heading}>
        { name }
      </Typography>

      <CardActions>
      <Grid direction='row' container>
        <Grid direction='column' container item xs={9}>
          <Grid direction='row' className={classes.statusline} container>
            <Grid container className={classes.title} item xs={4}>
              <Typography>Last Ran</Typography>
            </Grid>
            <Grid item xs={8}>
              <Typography>{ran_at ? formatTime(ran_at) : 'never'}</Typography>
            </Grid>
          </Grid>
          <Grid direction='row' className={classes.statusline} container>
            <Grid container className={classes.title} item xs={4}>
              <Typography>Next Run</Typography>
            </Grid>
            <Grid item xs={8}>
              <Typography>{runTime(ran_at, run_interval)} </Typography>
            </Grid>
          </Grid>
          <Grid direction='row' className={classes.statusline} container>
            <Grid container className={classes.title} item xs={4}>
              <Typography>Status</Typography>
            </Grid>
            <Grid className={classes[status]} direction='row' item container xs={8}>
              <Grid item><StatusIcon status={status}/></Grid>
              <Grid item><Typography>{ status }</Typography></Grid>
            </Grid>
          </Grid>
        </Grid>
        <Grid item container
          direction='row'
          className={classes.buttons}
          alignItems='center'
          xs={3}>
          <EtlButton mode='run' onClick={ toggleMode }><PlayArrowIcon/></EtlButton>
          <EtlButton mode='logs' onClick={ toggleMode }><LibraryBooksIcon/></EtlButton>
          <EtlButton mode='configure' onClick={ toggleMode }><SettingsIcon/></EtlButton>
          <EtlButton mode='remove' onClick={ toggleMode }><DeleteIcon/></EtlButton>
        </Grid>
      </Grid>
      </CardActions>
      <Collapse in={mode == 'configure'} timeout='auto' unmountOnExit>
        <CardContent>
          <ConfigScript script={config} schema={schema } />
        </CardContent>
      </Collapse>
      <Collapse in={mode == 'run'} timeout='auto' unmountOnExit>
        <CardContent>
          <SetRunState run_interval={run_interval} update={ update => postUpdate(update) } />
        </CardContent>
      </Collapse>
      <Collapse in={mode == 'remove'} timeout='auto' unmountOnExit>
        <CardContent>
          <Button onClick={ () => postUpdate({ archived: true }) }>Remove</Button>
        </CardContent>
      </Collapse>
      <Collapse in={mode == 'logs'} timeout='auto' unmountOnExit>
        <CardContent>
          <div className={classes.editor}>
            <Controlled
              options = {{
                readOnly: true,
                lineNumbers: true,
                lineWrapping: true,
                gutters: ['CodeMirror-lint-markers'],
                lint: true,
                tabSize: 2
              }}
              value={output}
              onBeforeChange={(editor, data, value) => { }}
            />
          </div>
        </CardContent>
      </Collapse>
    </CardContent>
  </Card>
}

export default EtlConfig;
