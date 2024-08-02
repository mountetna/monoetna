import React, { useContext } from 'react';
import { makeStyles } from '@mui/styles';

import AntSwitch from './ant_switch';
import {QueryResultsContext} from '../../contexts/query/query_results_context';
import Grid from '@mui/material/Grid';
import Button from '@mui/material/Button';
import Dialog from '@mui/material/Dialog';
import DialogActions from '@mui/material/DialogActions';
import DialogContent from '@mui/material/DialogContent';
import DialogContentText from '@mui/material/DialogContentText';
import DialogTitle from '@mui/material/DialogTitle';

const useStyles = makeStyles({
  table: {
    minWidth: 650
  },
  table_controls: {
    padding: '0px 15px'
  }
});

const QueryOptions = ({
  open, onClose
}: {
  open: boolean;
  onClose: Function;
}) => {
  const classes = useStyles();
  const {
    state: {
      expandMatrices,
      flattenQuery
    },
    setExpandMatrices,
    setFlattenQuery
  } = useContext(QueryResultsContext);

  return (
      <Dialog
        open={open}
        onClose={onClose}
        aria-labelledby="alert-dialog-title"
        aria-describedby="alert-dialog-description"
      >
        <DialogTitle id="alert-dialog-title">
          Query Options
        </DialogTitle>
        <DialogContent>
          <AntSwitch
            checked={expandMatrices}
            onChange={() => setExpandMatrices(!expandMatrices)}
            name='expand-matrices-query'
            leftOption='Nest matrices'
            rightOption='Expand matrices'
          />
          <AntSwitch
            checked={flattenQuery}
            onChange={() => setFlattenQuery(!flattenQuery)}
            name='flatten-query'
            leftOption='Nested'
            rightOption='Flattened'
          />
        </DialogContent>
      </Dialog>
  );
};

export default QueryOptions;
