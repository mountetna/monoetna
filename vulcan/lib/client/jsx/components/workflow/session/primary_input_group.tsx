import React, {useState, useContext, useEffect, useCallback, useMemo} from 'react';

import {VulcanContext} from '../../../contexts/vulcan_context';

import UserInput from '../user_interactions/inputs/user_input';
import {bindInputSpecification, InputSpecification} from '../user_interactions/inputs/input_types';

import AccordionSummary from '@material-ui/core/AccordionSummary';
import Grid from '@material-ui/core/Grid';
import AccordionDetails from '@material-ui/core/AccordionDetails';
import Accordion from '@material-ui/core/Accordion';
import Typography from '@material-ui/core/Typography';
import {makeStyles} from '@material-ui/core/styles';
import { BufferedInputsContext, WithBufferedInputs } from '../../../contexts/input_state_management';

const useStyles = makeStyles((theme) => ({
  header: {
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

function PrimaryInputInner({input}: {
  input: InputSpecification;
}) {
  const {state} = useContext(VulcanContext);
  const {status, workspace} = state;
  if (!workspace || !workspace.vulcan_config) return null;
  const {values, setValues} = useContext(BufferedInputsContext);

  const bound = useMemo(() => bindInputSpecification(
      input,
      // workspace.steps,
      workspace.vulcan_config,
      status.last_params,
      status.file_contents,
      status.params,
      status.ui_contents,
      values,
      setValues
    ), [
      input,
      values,
      setValues,
      // workspace.steps,
      workspace.vulcan_config,
      status.last_params,
      status.file_contents,
      status.params,
      workspace
    ])

  return <UserInput
    input={bound}
  />
}

export default function PrimaryInputGroup({inputs, groupName, select, expanded}: {
  inputs: InputSpecification[],
  groupName: string,
  expanded: boolean,
  select: (event: any) => void
}) {
  const {commitSessionInputChanges, dispatch} = useContext(VulcanContext);
  const sortedInputs = useMemo(() => inputs.sort(
    (a, b) =>
      collator.compare(a.label, b.label)), [inputs]);

  const classes = useStyles();

  return (<Accordion elevation={0} expanded={ expanded } onChange={ select }>
    <AccordionSummary className={classes.header}>
        <Typography>{groupName}</Typography>
    </AccordionSummary>
    <AccordionDetails className={classes.group}>
      <Grid container direction='column'>
      {sortedInputs.map((input, index) => {
        return (
          <WithBufferedInputs
            commitSessionInputChanges={commitSessionInputChanges}
            dispatch={dispatch}
            stepName={input.name}
            key={index}
          >
            <PrimaryInputInner
              input={input}
            />
          </WithBufferedInputs>
        );
      })}
      </Grid>
    </AccordionDetails>
  </Accordion>);
}
