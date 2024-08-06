import React, { useContext } from 'react';
import { makeStyles } from '@mui/styles';

import {QueryResultsContext} from '../../contexts/query/query_results_context';
import Grid from '@mui/material/Grid';
import FormControlLabel from '@mui/material/FormControlLabel';
import FormGroup from '@mui/material/FormGroup';
import Checkbox from '@mui/material/Checkbox';
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
      flattenQuery,
      showDisconnected
    },
    setExpandMatrices,
    setFlattenQuery,
    setShowDisconnected
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
          <FormGroup>
            <FormControlLabel
              label="Expand matrix attributes to columns"
              control={
                <Checkbox
                  checked={expandMatrices}
                  onChange={() => setExpandMatrices(!expandMatrices)}
                />
              }
            />
            <FormControlLabel
              label="Flatten nested attributes"
              control={
                <Checkbox
                  checked={flattenQuery}
                  onChange={() => setFlattenQuery(!flattenQuery)}
                />
              }
            />
            <FormControlLabel
              label="Show only rows from disconnected records"
              control={
                <Checkbox
                  checked={showDisconnected}
                  onChange={() => setShowDisconnected(!showDisconnected)}
                />
              }
            />
          </FormGroup>
        </DialogContent>
      </Dialog>
  );
};

export default QueryOptions;
