import React, {useState, useCallback} from 'react';

import TableRow from '@material-ui/core/TableRow';
import TableCell from '@material-ui/core/TableCell';
import {makeStyles} from '@material-ui/core/styles';

import {formatTime, runTime, runDesc} from './run-state';

import {Status, Job} from '../polyphemus';

const useStyles = makeStyles((theme) => ({
  workflowrow: {
    marginBottom: '15px',
    cursor: 'pointer',
    '&:hover': {
      background: '#eee'
    }
  }
}));

const WorkflowStatus = ({
  status,
  onClick
}:{ status: Status; onClick: Function}) => {
  const { workflow_name, workflow_type, pipeline_state, pipeline_finished_at } = status;

  const classes: any = useStyles();

  return (
    <TableRow className={classes.workflowrow} onClick={ onClick as React.MouseEventHandler<HTMLTableRowElement>}>
      <TableCell>{workflow_type}</TableCell>
      <TableCell>{workflow_name}</TableCell>
      <TableCell align="right">{pipeline_state || 'never run'}</TableCell>
      <TableCell align="right">{pipeline_finished_at ? formatTime(pipeline_finished_at) : 'never'}</TableCell>
    </TableRow>
  );
};

export default WorkflowStatus
