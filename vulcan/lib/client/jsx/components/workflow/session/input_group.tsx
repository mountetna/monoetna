import React, {useState, useContext, useEffect, useCallback, useMemo} from 'react';

import {VulcanContext} from '../../../contexts/vulcan_context';

import UserInput from '../user_interactions/inputs/user_input';
import {
  completedStepNames, paramUINames, pendingStepNames, pendingUIInputStepReady, uiComponentOfStep
} from '../../../selectors/workflow_selectors';
import {BoundInputSpecification} from '../user_interactions/inputs/input_types';
import {useWorkspace} from '../../../contexts/workspace_context';

import IconButton from '@material-ui/core/IconButton';
import AccordionSummary from '@material-ui/core/AccordionSummary';
import Grid from '@material-ui/core/Grid';
import AccordionDetails from '@material-ui/core/AccordionDetails';
import Accordion from '@material-ui/core/Accordion';
import Typography from '@material-ui/core/Typography';
import {makeStyles} from '@material-ui/core/styles';

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

interface Props {
  inputs: BoundInputSpecification[],
  groupName: string,
  expanded: boolean,
  select: (event: any) => void
}

const collator = new Intl.Collator(undefined, {
  numeric: true, sensitivity: 'base'
});

export default function InputGroup({inputs, groupName, select, expanded}: Props) {

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
        return (<UserInput
          input={input}
          key={index}
        />);
      })}
      </Grid>
    </AccordionDetails>
  </Accordion>);
}
