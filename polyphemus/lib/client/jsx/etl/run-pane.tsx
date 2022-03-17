import React, {useState, useEffect, useCallback, useMemo} from 'react';
import * as _ from 'lodash';

import Button from '@material-ui/core/Button';
import Grid from '@material-ui/core/Grid';
import MenuItem from '@material-ui/core/MenuItem';
import Select from '@material-ui/core/Select';
import Switch from '@material-ui/core/Switch';
import TextField from '@material-ui/core/TextField';
import Typography from '@material-ui/core/Typography';
import InputAdornment from '@material-ui/core/InputAdornment';
import {makeStyles} from '@material-ui/core/styles';
import EtlPane, {EtlPaneHeader} from './etl-pane';

import {
  RUN_ONCE,
  RUN_NEVER,
  RUN_INTERVAL,
  getRunState,
  getRunIntervalTime,
  runTime
} from './run-state';

const useStyles = makeStyles((theme) => ({
  runbar: {
    width: 'auto',
    justifyContent: 'space-between'
  },
  title: {
    flex: '0 0 200px'
  },
  values: {
    flex: '1 1 auto'
  },
  params: {
    flex: '1 1 auto'
  },
  param: {
    borderBottom: '1px solid #eee'
  }
}));

type Value = undefined | string | boolean;

type SelectParam = {
  value: string;
  label: string;
};

const Param = ({
  value,
  opts,
  name,
  update
}: {
  value: Value;
  opts: SelectParam[] | 'string' | 'boolean';
  name: string;
  update: (name: string, value: string | boolean) => void;
}) => {
  if (Array.isArray(opts)) {
    return (
      <Select
        value={value || ''}
        onChange={(e) => update(name, e.target.value as string)}
      >
        {opts.map((opt) => (
          <MenuItem key={opt.value} value={opt.value}>
            {`${opt.value}${opt.label ? ` - ${opt.label}` : ''}`}
          </MenuItem>
        ))}
      </Select>
    );
  }

  if (opts == 'string') {
    return (
      <TextField
        size='small'
        fullWidth
        value={value == undefined ? '' : value}
        onChange={(e) => update(name, e.target.value)}
      />
    );
  }

  if (opts == 'boolean') {
    return (
      <Switch
        size='small'
        checked={value == undefined ? false : (value as boolean)}
        onChange={(e) => update(name, e.target.checked)}
      />
    );
  }

  return null;
};

const RunPane = ({
  run_interval,
  update,
  params,
  param_opts,
  selected,
  error
}: {
  run_interval: number;
  params: any;
  param_opts: any;
  update: Function;
  error: string;
  selected: string | null;
}) => {
  const [runState, setRunState] = useState(getRunState(RUN_ONCE));
  const [runIntervalTime, setRunIntervalTime] = useState(
    getRunIntervalTime(RUN_ONCE)
  );
  const [newParams, setParams] = useState<any>({});

  useEffect(() => {
    setParams(params);
  }, [params]);

  const runValue = () =>
    runState == RUN_INTERVAL ? runIntervalTime : runState;
  const reset = () => {
    setRunState(getRunState(RUN_ONCE));
    setRunIntervalTime(getRunIntervalTime(RUN_ONCE));
    setParams(params);
  };

  const classes = useStyles();

  const nonBooleanParamOpts = useMemo(() => {
    return Object.entries(param_opts)
      .filter(([param, value]) => value !== 'boolean')
      .map(([param, value]) => param)
      .sort();
  }, [param_opts]);

  const formValid = useMemo(() => {
    return (
      (_.isEqual(nonBooleanParamOpts, Object.keys(newParams).sort()) ||
        _.isEqual(
          Object.keys(param_opts).sort(),
          Object.keys(newParams).sort()
        )) &&
      Object.values(newParams).every((v) => v != null && v !== '')
    );
  }, [nonBooleanParamOpts, param_opts, newParams]);

  return (
    <EtlPane mode='run' selected={selected}>
      <EtlPaneHeader title='Run' error={error}>
        <Grid className={classes.runbar} spacing={1} container>
          <Grid item>
            <Button
              onClick={() =>
                update({run_interval: runValue(), params: newParams})
              }
              disabled={!formValid}
            >
              Save
            </Button>
          </Grid>
          <Grid item>
            <Button onClick={reset} color='secondary'>
              Reset
            </Button>
          </Grid>
        </Grid>
      </EtlPaneHeader>
      <Grid spacing={1} container alignItems='center'>
        <Grid className={classes.title} item>
          <Typography>Interval</Typography>
        </Grid>
        <Grid className={classes.values} item>
          <Select
            value={runState}
            onChange={(e) => setRunState(parseInt(e.target.value as string))}
          >
            <MenuItem value={RUN_ONCE}>Once</MenuItem>
            <MenuItem value={RUN_NEVER}>Never</MenuItem>
            <MenuItem value={RUN_INTERVAL}>Every</MenuItem>
          </Select>
        </Grid>
        {runState == RUN_INTERVAL && (
          <Grid item>
            {' '}
            <TextField
              size='small'
              InputProps={{
                endAdornment: (
                  <InputAdornment position='end'>seconds</InputAdornment>
                )
              }}
              value={runIntervalTime}
              onChange={(e) => {
                setRunIntervalTime(parseInt(e.target.value as string));
              }}
            />
          </Grid>
        )}
      </Grid>
      <Grid spacing={1} container alignItems='flex-start'>
        <Grid className={classes.title} item>
          <Typography>Parameters</Typography>
        </Grid>
        <Grid className={classes.params} item>
          {Object.keys(param_opts)
            .sort()
            .map((param_name) => (
              <Grid
                alignItems='center'
                key={param_name}
                className={classes.param}
                container
              >
                <Grid item xs={4}>
                  {param_name}
                </Grid>
                <Grid item xs={5}>
                  <Param
                    name={param_name}
                    value={newParams && (newParams[param_name] as Value)}
                    opts={param_opts[param_name]}
                    update={(name, value) => {
                      console.log({name, value});
                      setParams({...newParams, [name]: value});
                    }}
                  />
                </Grid>
              </Grid>
            ))}
        </Grid>
      </Grid>
    </EtlPane>
  );
};

export default RunPane;
