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
import ReportProblemOutlinedIcon from '@material-ui/icons/ReportProblemOutlined';

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

type Value = undefined | string | boolean | number;

type SelectValue = {
  value: string;
  description: string;
};

type ParamType = 'select' | 'options' | 'boolean' | 'integer' | 'string' | 'datetime';

type ComplexParam = {
  type: ParamType;
  values?: SelectValue[];
  valuesFrom?: string;
  description?: string;
  default?: string;
};

type ParamComponent = {
  value: Value;
  options: ComplexParam;
  name: string;
  update: (name: string, value: Value) => void;
  config: any;
};

const SelectParam = ({value, options,name,update,config}:ParamComponent) => {
  const { values, default: defaultValue } = options;

  if (!values) return null;

  return (
    <Select
      value={value ? value : defaultValue ? defaultValue : ''}
      onChange={(e) => update(name, e.target.value as string)}
    >
      {
        values.map( opt =>
          <MenuItem key={opt.value} value={opt.value}>
            {`${opt.value}${opt.description ? ` - ${opt.description}` : ''}`}
          </MenuItem>
        )
      }
    </Select>
  );
};

const OptionsParam = ({value, options,name,update,config}:ParamComponent) => {
  const { values, valuesFrom } = options;
  let onChange: (e: any, v: string[]) => void;
  let fillValues: string[];

  if (valuesFrom === 'model_names') {
    fillValues = ['all'].concat(Object.keys(config).sort());
    onChange = (e, v) => update(name, v.includes('all') ? 'all' : v.join(','));
  } else if (values) {
    fillValues = values.map(v => v.value);
    onChange = (e, v) => update(name, v.join(','));
  } else {
    return null;
  }

  return (
    <Autocomplete
      multiple
      fullWidth
      id={name}
      options={fillValues}
      value={'' === value ? [] : (value as string)?.split(',') || []}
      renderInput={(params) => <TextField {...params} variant='standard' />}
      onChange={onChange}
    />
  );
}

const StringParam = ({value, options,name,update,config}:ParamComponent) => {
  return (
    <TextField
      size='small'
      placeholder={options.description}
      fullWidth
      value={value == undefined ? '' : value}
      onChange={(e) => update(name, e.target.value)}
    />
  );
}

const BooleanParam = ({value, options,name,update,config}:ParamComponent) => {
  return (
    <Switch
      size='small'
      checked={value == undefined ? false : (value as boolean)}
      onChange={(e) => update(name, e.target.checked)}
    />
  );
}

const IntegerParam = ({value, options,name,update,config}:ParamComponent) => {
  return (
    <TextField
      size='small'
      type='number'
      placeholder={options.description}
      fullWidth
      value={value == undefined ? '' : value}
      onChange={(e) => update(name, parseInt(e.target.value))}
    />
  );
}

const DateTimeParam = ({value, options,name,update,config}:ParamComponent) => {
  return (
    <TextField
      size='small'
      type='datetime-local'
      placeholder={options.description}
      fullWidth
      value={value == undefined ? '' : value}
      onChange={(e) => { console.log({ v: e.target.value}); update(name, e.target.value)} }
    />
  );
}

const PARAMS = {
  select: SelectParam,
  options: OptionsParam,
  boolean: BooleanParam,
  integer: IntegerParam,
  string: StringParam,
  datetime: DateTimeParam,
}

const Param = ({
  value,
  options,
  name,
  update,
  config
}: ParamComponent) => {
  if (options.type in PARAMS) {
    const Component = PARAMS[options.type];

    return <Component
      value={value}
      options={options}
      name={name}
      update={update}
      config={config}
    />;
  }

  return null;
};

const paramsWithDefaults = (params: any, param_opts: { [param_name: string]: ComplexParam }) => {
  let newParams = {...params};
  Object.entries(param_opts).forEach(
    ([param_name, opts]) => {
      if (param_name in newParams) return;

      if ('default' in opts) newParams[param_name] = opts.default;
      else if (opts.type == 'boolean') newParams[param_name] = false;
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

  const formValid = useMemo(() => {
    const newParamOpts = Object.keys(paramsWithDefaults(newParams, param_opts)).sort();
    const oldParamOpts = Object.keys(param_opts).sort();

    return _.isEqual(oldParamOpts, newParamOpts);
  }, [param_opts, newParams]);

  const savePayload = useCallback(() => {
    return {run_interval: runIntervalTime, config: paramsWithDefaults(newParams, param_opts) };
  }, [runIntervalTime, newParams]);

  const formChanged = useMemo(() => {
    return !_.isEqual(
      savePayload(),
      {
        config: paramsWithDefaults(params, param_opts),
        run_interval: getRunIntervalTime(run_interval)
      }
    );
  }, [savePayload, params, param_opts]);

  const dryRun = useMemo(() => {
    const nextParams = {...params, ...newParams}
    return 'commit' in nextParams && !nextParams['commit']
  }, [param_opts, params, newParams])

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
          { dryRun && 
            <Grid item>
              <Tooltip title="Parsed updates will not be committed to the database">
                <Grid item container alignItems='center' spacing={0}>
                  <Grid item>
                    <ReportProblemOutlinedIcon/>
                  </Grid>
                  <Grid item>
                    <Typography>Dry Run Mode</Typography>
                  </Grid>
                </Grid>
              </Tooltip>
            </Grid>
          }
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
                    options={param_opts[param_name]}
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
