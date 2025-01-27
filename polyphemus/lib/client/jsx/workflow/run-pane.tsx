import React, {useState, useEffect, useCallback, useMemo} from 'react';
import * as _ from 'lodash';
import {json_post} from 'etna-js/utils/fetch';

import Button from '@material-ui/core/Button';
import Tooltip from '@material-ui/core/Tooltip';
import Grid from '@material-ui/core/Grid';
import MenuItem from '@material-ui/core/MenuItem';
import Select from '@material-ui/core/Select';
import Switch from '@material-ui/core/Switch';
import TextField from '@material-ui/core/TextField';
import Typography from '@material-ui/core/Typography';
import InputAdornment from '@material-ui/core/InputAdornment';
import {makeStyles} from '@material-ui/core/styles';
import WorkflowPane, {WorkflowPaneHeader} from './workflow-pane';
import Autocomplete from '@material-ui/lab/Autocomplete';

import {
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
  completed: {
    color: 'green',
    paddingLeft: '0.5rem'
  },
  error: {
    paddingLeft: '0.5rem',
    color: 'red'
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
  description: string;
  default?: boolean;
};

type ComplexParam = {
  type: string;
  value: string;
};

const DefaultSelect = ({
  value,
  opts,
  name,
  update
}: {
  value: Value;
  opts: SelectParam[];
  name: string;
  update: (name: string, value: string | boolean) => void;
}) => {
  const defaultOption = opts.find((o) => o.default);

  return (
    <Select
      value={value ? value : defaultOption ? defaultOption.value : ''}
      onChange={(e) => update(name, e.target.value as string)}
    >
      {opts.map((opt) => (
        <MenuItem key={opt.value} value={opt.value}>
          {`${opt.value}${opt.description ? ` - ${opt.description}` : ''}`}
        </MenuItem>
      ))}
    </Select>
  );
};

const Param = ({
  value,
  opts,
  name,
  update,
  config
}: {
  value: Value;
  opts: SelectParam[] | ComplexParam | 'string' | 'boolean';
  name: string;
  update: (name: string, value: string | boolean) => void;
  config: any;
}) => {
  if (Array.isArray(opts)) {
    return (
      <DefaultSelect value={value} opts={opts} name={name} update={update} />
    );
  } else if (_.isObject(opts)) {
    if (opts.type && opts.type === 'options' && opts.value === 'model_names') {
      const modelNames = ['all'].concat(Object.keys(config).sort());
      return (
        <Autocomplete
          multiple
          fullWidth
          id={name}
          options={modelNames}
          value={'' === value ? [] : (value as string)?.split(',') || []}
          renderInput={(params) => <TextField {...params} variant='standard' />}
          onChange={(e, v) => {
            if (v.includes('all')) {
              update(name, 'all');
            } else {
              update(name, v.join(','));
            }
          }}
        />
      );
    } else {
      return (
        <TextField
          size='small'
          fullWidth
          value={value == undefined ? '' : value}
          onChange={(e) => update(name, e.target.value)}
        />
      );
    }
  } else if (opts == 'string') {
    return (
      <TextField
        size='small'
        fullWidth
        value={value == undefined ? '' : value}
        onChange={(e) => update(name, e.target.value)}
      />
    );
  } else if (opts == 'boolean') {
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

const paramsWithDefaults = (params, param_opts) => {
  let newParams = {...params};
  Object.entries(param_opts).forEach(
    ([param_name,opts]) => {
      if (param_name in newParams) return;
      if (!Array.isArray(opts)) return;

      const defaultOption = opts.find((o) => o.default);

      if (defaultOption) {
        newParams[param_name] = defaultOption.value;
      }
    }
  );
  return newParams;
}

const RunPane = ({
  project_name,
  config_id,
  run_interval,
  update,
  params,
  param_opts,
  selected,
  config
}: {
  project_name: string;
  config_id: number;
  run_interval: number;
  params: any;
  param_opts: any;
  update: Function;
  selected: string | null;
  config: any;
}) => {
  const [runIntervalTime, setRunIntervalTime] = useState(getRunIntervalTime(run_interval));
  const [newParams, setParams] = useState<any>({});
  const [error, setError] = useState('');
  const [message, setMessage] = useState('');
  const [ranOnce, setRanOnce] = useState(false)

  useEffect(() => setParams(paramsWithDefaults(params, param_opts)), [params, param_opts]);
  useEffect(() => setRunIntervalTime(getRunIntervalTime(run_interval)), [run_interval]);

  const reset = useCallback(() => {
      setError('');
      setMessage('');
      setRunIntervalTime(getRunIntervalTime(run_interval));
      setParams(paramsWithDefaults(params, param_opts));
    }, [params, param_opts]
  );

  useEffect( () => reset(), []);

  const classes = useStyles();

  const nonBooleanParamOpts = useMemo(() => {
    return Object.entries(param_opts)
      .filter(([param, value]) => value !== 'boolean')
      .map(([param, value]) => param)
      .sort();
  }, [param_opts]);

  const formValid = useMemo(() => {
    const newParamOpts = Object.keys(newParams).sort();
    const nonBoolFilled = _.isEqual(nonBooleanParamOpts, newParamOpts);
    const allFilled = _.isEqual(Object.keys(param_opts).sort(), newParamOpts);
    const noneEmpty = Object.values(newParams).every(
      (v) => v != null && v !== '' && !_.isEqual(v, [])
    );

    return (nonBoolFilled || allFilled) && noneEmpty;
  }, [nonBooleanParamOpts, param_opts, newParams]);

  const savePayload = useCallback(() => {
    return {run_interval: runIntervalTime, config: newParams};
  }, [runIntervalTime, newParams]);

  const formChanged = useMemo(() => {
    console.log({ sp: savePayload(), np: { config: paramsWithDefaults(params, param_opts), run_interval } });
    return !_.isEqual(savePayload(), { config: paramsWithDefaults(params, param_opts), run_interval: getRunIntervalTime(run_interval) });
  }, [savePayload, params, param_opts]);

  console.log({nonBooleanParamOpts, newParams, param_opts});
  console.log({formChanged, formValid});

  return (
    <WorkflowPane mode='run' selected={selected}>
      <WorkflowPaneHeader title='Run'>
        <Grid className={classes.runbar} spacing={1} container alignItems='center'>
          <Grid item>
            <Tooltip title={ !formValid ? 'Missing values' : !formChanged ? 'No changes' : '' }>
              <span>
                <Button
                  disabled={ !formChanged || !formValid }
                  onClick={() => {
                    if (runIntervalTime && runIntervalTime < 300) {
                      setError('Run interval cannot be less than 300 seconds.')
                    } else {
                      setError('');
                      update(savePayload());
                    }
                  }}
                >
                  Save
                </Button>
              </span>
            </Tooltip>
          </Grid>
          <Grid item>
            <Tooltip title={ !formChanged ? 'No changes' : '' }>
              <span>
                <Button onClick={reset} color='secondary' disabled={ !formChanged }>
                  Reset
                </Button>
              </span>
            </Tooltip>
          </Grid>
          <Grid item>
            <Tooltip title={ ranOnce ? 'Already ran' : !formValid ? 'Missing values' : '' }>
              <span>
                <Button
                  disabled={ !formValid || ranOnce}
                  onClick={() => {
                    setRanOnce(true);
                    json_post(`/api/workflows/${project_name}/runtime_configs/run_once/${config_id}`, { config: newParams })
                      .then((run) => setMessage('Launched!'))
                      .catch((r) => r.then(({error}: {error: string}) => setError(error)));
                  }}
                >
                  Run Once
                </Button>
              </span>
            </Tooltip>
          </Grid>
          {(error || message) && (
            <Grid item className={error ? classes.error : classes.completed}>
              <Typography>{error || message}</Typography>
            </Grid>
          )}
        </Grid>
      </WorkflowPaneHeader>
      <Grid spacing={1} container alignItems='center'>
        <Grid className={classes.title} item>
          <Typography>Interval</Typography>
        </Grid>
        <Grid item xs={4}>
          {' '}
          <TextField
            size='small'
            InputProps={{
              endAdornment: (
                <InputAdornment position='end'>{'seconds (300 minimum, 0 to disable)'}</InputAdornment>
              )
            }}
            value={runIntervalTime}
            onChange={(e) => {
              setRunIntervalTime(parseInt(e.target.value as string || '0'));
            }}
          />
        </Grid>
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
                    config={config}
                    name={param_name}
                    value={newParams && (newParams[param_name] as Value)}
                    opts={param_opts[param_name]}
                    update={(name, value) => {
                      setParams({...newParams, [name]: value});
                    }}
                  />
                </Grid>
              </Grid>
            ))}
        </Grid>
      </Grid>
    </WorkflowPane>
  );
};

export default RunPane;
