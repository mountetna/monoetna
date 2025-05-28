import React, { useState, useContext, useCallback } from 'react';
import { makeStyles } from '@material-ui/core/styles';

import {QueryResultsContext} from '../../contexts/query/query_results_context';
import {SavedQuery} from '../../contexts/query/query_types';
import Grid from '@material-ui/core/Grid';
import FormControlLabel from '@material-ui/core/FormControlLabel';
import FormGroup from '@material-ui/core/FormGroup';
import Checkbox from '@material-ui/core/Checkbox';
import Button from '@material-ui/core/Button';
import Dialog from '@material-ui/core/Dialog';
import DialogActions from '@material-ui/core/DialogActions';
import DialogContent from '@material-ui/core/DialogContent';
import DialogTitle from '@material-ui/core/DialogTitle';
import TextField from '@material-ui/core/TextField';
import {json_post} from 'etna-js/utils/fetch';

const useStyles = makeStyles({
  table: {
    minWidth: 650
  },
  table_controls: {
    padding: '0px 15px'
  }
});

const QuerySaveModal = ({
  open, onClose
}: {
  open: boolean;
  onClose: () => void;
}) => {
  const classes = useStyles();
  const { state: { savedQueries }, setSavedQueries } = useContext(QueryResultsContext);

  const [ comment, setComment ] = useState('');

  const handleClose = useCallback(
    () => {
      setComment('');
      onClose();
    }, [onClose]
  );

  const handleSave = useCallback(
    () => {
      const params = new URLSearchParams(window.location.hash.slice(1));

      if (!params.has('q')) {
        handleClose();
        return;
      }

      const query = params.get('q');

      json_post(`/api/query_history/${CONFIG.project_name}/create`, {
        comment,
        query
      }).then( ({query}:{query: SavedQuery}) => {
        setSavedQueries(savedQueries.concat(query));
      });

      handleClose();
    }, [comment]
  );

  return (
      <Dialog
        open={open}
        onClose={onClose}
        fullWidth
        maxWidth='sm'
        aria-labelledby="alert-dialog-title"
        aria-describedby="alert-dialog-description"
      >
        <DialogTitle id="alert-dialog-title">
          Save Query
        </DialogTitle>
        <DialogContent>
          <TextField
            fullWidth
            value={comment}
            placeholder='Query description'
            onChange={(e: React.ChangeEvent<any>) => setComment(e.target.value)}
          />
        </DialogContent>
        <DialogActions>
          <Button color='secondary' onClick={handleClose}>Cancel</Button>
          <Button color='primary' disabled={!comment} onClick={handleSave}>Save</Button>
        </DialogActions>
      </Dialog>
  );
};

export default QuerySaveModal;
