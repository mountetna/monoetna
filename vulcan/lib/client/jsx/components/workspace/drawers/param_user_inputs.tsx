import React, {
  useState,
  useContext,
  useMemo
} from 'react';

import {VulcanContext} from '../../../contexts/vulcan_context';

import {
  inputGroupName,
} from '../../../selectors/workflow_selectors';
import {
  getParamUISpecifications,
  bindInputSpecification,
  InputSpecification
} from '../ui_definitions/input_types';
import { BufferedInputsContext, WithBufferedInputs } from '../../../contexts/input_state_management';
import UserInput from './user_input';

import Collapse from '@material-ui/core/Collapse';
import Typography from '@material-ui/core/Typography';
import Grid from '@material-ui/core/Grid';
import Card from '@material-ui/core/Card';
import ExpandMoreIcon from '@material-ui/icons/ExpandMore';
import ExpandLessIcon from '@material-ui/icons/ExpandLess';
import IconButton from '@material-ui/core/IconButton';
import AccordionSummary from '@material-ui/core/AccordionSummary';
import AccordionDetails from '@material-ui/core/AccordionDetails';
import Accordion from '@material-ui/core/Accordion';
import {makeStyles} from '@material-ui/core/styles';
import Tooltip from '@material-ui/core/Tooltip';

const useStyles = makeStyles((theme) => ({
  card: {
    borderRadius: 0,
    border: '1px solid #eee'
  },
  header: {
    padding: '5px 10px',
    cursor: 'pointer'
  },
  groupHeader: {
    background: '#eee',
    cursor: 'pointer',
    minHeight: '32px',
    height: '32px'
  },
  group: {
    padding: 0,
    borderBottom: '1px solid #eee'
  }
}));

const collator = new Intl.Collator(undefined, {
  numeric: true, sensitivity: 'base'
});

function ParamInputIndividual({input}: {
  input: InputSpecification;
}) {
  const {state} = useContext(VulcanContext);
  const {status, workspace} = state;
  if (!workspace || !workspace.vulcan_config) return null;
  const {values, setValues, showError} = useContext(BufferedInputsContext);

  const bound = useMemo(() => bindInputSpecification(
      input,
      workspace.vulcan_config,
      status.last_params,
      status.file_contents,
      status.params,
      status.ui_contents,
      values,
      setValues,
      showError
    ), [
      input,
      values,
      setValues,
      workspace.vulcan_config,
      status.last_params,
      status.file_contents,
      status.params,
      workspace
    ])

  return <UserInput
    input={bound} hideLabel={true}
  />
}

function ParamInputGroup({inputs, groupName, select, expanded}: {
  inputs: InputSpecification[],
  groupName: string,
  expanded: boolean,
  select: (event: any) => void
}) {
  const {commitSessionInputChanges, dispatch, useActionInvoker} = useContext(VulcanContext);
  const sortedInputs = useMemo(() => inputs.sort(
    (a, b) =>
      collator.compare(a.label, b.label)), [inputs]);

  const classes = useStyles();

  return (<Accordion elevation={0} expanded={ expanded } onChange={ select }>
    <AccordionSummary className={classes.groupHeader}>
        <Typography>{groupName}</Typography>
    </AccordionSummary>
    <AccordionDetails className={classes.group}>
      <Grid container direction='column'>
      {sortedInputs.map((input, index) => {
        return (
          <WithBufferedInputs
            commitSessionInputChanges={commitSessionInputChanges}
            dispatch={dispatch}
            invoke={useActionInvoker()}
            stepName={input.name}
            key={index}
          >
            <ParamInputIndividual
              input={input}
            />
          </WithBufferedInputs>
        );
      })}
      </Grid>
    </AccordionDetails>
  </Accordion>);
}

export default function ParamInputs() {
  const {state} = useContext(VulcanContext);
  const {workspace} = state;
  if (!workspace) return null;

  const previousParams = Object.keys(state.status.last_params).length==0;

  const inputSpecifications = useMemo(
    () => getParamUISpecifications(workspace),
    [workspace]
  );

  let groupedInputs = useMemo(() => {
    return inputSpecifications.reduce((result, spec) => {
      let groupName = inputGroupName(spec.label);
      if (groupName != null) {spec.label = spec.label.replace(`${groupName}__`, '')}
      groupName = groupName || 'Ungrouped';
      result[groupName] = result[groupName] || [];
      result[groupName].push(spec);
      return result;
    }, {} as {[k: string]: InputSpecification[]});
  }, [
    inputSpecifications
  ]);

  const [expanded, setExpanded] = useState(previousParams || Object.keys(groupedInputs).length==1 ? Object.keys(groupedInputs) : [] as string[]);
  function toggleExpanded(groupName: string) {
    if (expanded.includes(groupName)) {
      setExpanded(expanded.filter(val => val != groupName));
    } else {
      setExpanded(expanded.concat(groupName))
    }
  }

  const [open, setOpen] = useState(previousParams);

  const classes = useStyles();

  return (
    <Card className={classes.card}>
      <Tooltip title={`${open ? 'hide' : 'show'} output-independent parameters`}>
        <Grid
          className={classes.header}
          container
          alignItems='center'
          justifyContent='space-between'
          onClick={() => setOpen(!open)}
        >
          <Grid item>
            <Typography variant='h6'>Analysis Setup</Typography>
          </Grid>
          <Grid item>
            <IconButton size='small'>
              {open ? (
                <ExpandLessIcon fontSize='small' />
              ) : (
                <ExpandMoreIcon fontSize='small' />
              )}
            </IconButton>
          </Grid>
        </Grid>
      </Tooltip>
      <Collapse in={open}>
        {Object.keys(groupedInputs)
          .sort()
          .map((groupName, index) => {
            return (
              <ParamInputGroup
                expanded={expanded.includes(groupName)}
                select={() => toggleExpanded(groupName)}
                groupName={groupName.replace(/^\d__/,'').replace('_',' ')}
                key={index}
                inputs={groupedInputs[groupName]}
              />
            );
          })}
      </Collapse>
    </Card>
  );
}
