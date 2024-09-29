import React, { useContext } from 'react';
import { makeStyles } from '@material-ui/core/styles';

import {QueryResultsContext} from '../../contexts/query/query_results_context';
import Grid from '@material-ui/core/Grid';
import FormControlLabel from '@material-ui/core/FormControlLabel';
import FormGroup from '@material-ui/core/FormGroup';
import Checkbox from '@material-ui/core/Checkbox';
import Button from '@material-ui/core/Button';
import Dialog from '@material-ui/core/Dialog';
import DialogActions from '@material-ui/core/DialogActions';
import DialogContent from '@material-ui/core/DialogContent';
import DialogContentText from '@material-ui/core/DialogContentText';
import DialogTitle from '@material-ui/core/DialogTitle';

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
