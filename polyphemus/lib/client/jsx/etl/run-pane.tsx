import React, { useState, useCallback } from 'react';
import Button from '@material-ui/core/Button';
import Grid from '@material-ui/core/Grid';
import MenuItem from '@material-ui/core/MenuItem';
import Select from '@material-ui/core/Select';
import TextField from '@material-ui/core/TextField';
import InputAdornment from '@material-ui/core/InputAdornment';
import {makeStyles} from '@material-ui/core/styles';
import EtlPane, {EtlPaneHeader} from './etl-pane';

import { RUN_ONCE, RUN_NEVER, RUN_INTERVAL, getRunState, getRunIntervalTime, runTime } from './run-state';

const useStyles = makeStyles( theme => ({
  runbar: {
    width: 'auto',
    justifyContent: 'space-between'
  }
}));

const RunPane = ({run_interval,update,selected}:{run_interval:number, update:Function,selected:string}) => {
  const [ runState, setRunState ] = useState(getRunState(run_interval));
  const [ runIntervalTime, setRunIntervalTime ] = useState(getRunIntervalTime(run_interval));

  const runValue = () => runState == RUN_INTERVAL ? runIntervalTime : runState;
  const reset = () => { setRunState(getRunState(run_interval)); setRunIntervalTime(getRunIntervalTime(run_interval)); };

  const classes = useStyles();

  return <EtlPane mode='run' selected={selected}>
    <EtlPaneHeader title='Run'>
      <Grid className={classes.runbar} spacing={1} container alignItems='center'>
        <Grid item><Select value={runState} onChange={e => setRunState(parseInt(e.target.value as string))}>
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
    </EtlPaneHeader>
  </EtlPane>
}

export default RunPane;
