import React, {useState, useEffect, useMemo, useContext} from 'react';

import { makeStyles } from '@material-ui/core/styles';
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
import Typography from '@material-ui/core/Typography';
import Collapse from '@material-ui/core/Collapse';

import {json_get, json_delete} from 'etna-js/utils/fetch';
import Editor from 'etna-js/components/editor';

import {QueryResultsContext} from '../../contexts/query/query_results_context';
import {QueryResponse, SavedQuery} from '../../contexts/query/query_types';

const useStyles = makeStyles((theme) => ({
  query_header: {
    borderBottom: '1px solid #ccc',
    padding: '5px 0px'
  },
  diff: {
    height: '200px',
    border: '1px solid #ccc',
    '& .cm-editor': {
      height: '100%',
      width: '100%'
    },
    '& div': {
      height: '100%',
      width: '100%'
    },
    '& .cm-editor div': {
      height: 'auto',
      width: 'auto'
    }
  },
  actions: {
    paddingRight: '25px'
  }
}));

const formatTime = (time: string) =>
  Intl.DateTimeFormat('en-US', {
    year: 'numeric',
    month: 'short',
    day: 'numeric',
    hour: 'numeric',
    minute: 'numeric'
  }).format(new Date(time));

interface Item {
  comment: string;
  [key: string]: string;
};

const QueryHistory = ({
  open, onClose
}: {
  open: boolean;
  onClose: () => void;
}) => {
  const classes = useStyles();

  const [selectedIndex, setIndex] = useState(-1);

  const {state: { savedQueries }, setSavedQueries, setQueryStateFromString} = useContext(QueryResultsContext);

  useEffect(
    () => {
      json_get(`/api/query_history/${CONFIG.project_name}`).then(
        ({queries}) => setSavedQueries(queries)
      )
    }, []
  );

  const removeQuery = (query: SavedQuery) => {
    json_delete(`/api/query_history/${CONFIG.project_name}/${query.id}`).then(
      () => {
        setSavedQueries(
          savedQueries.filter( (q,i) => i != selectedIndex )
        );
        setIndex(-1);
      }
    );
  };

  const handleClose = () => {
    setIndex(-1);
    onClose();
  }

  const loadQuery = (query: SavedQuery) => {
    const { pathname } = window.location;
    setQueryStateFromString( `q=${query.query}`);
    handleClose();
  };

  return (
    <Dialog
      maxWidth='md'
      fullWidth
      scroll='paper'
      open={open}
      onClose={handleClose}
      aria-labelledby='form-dialog-title'
    >
      <DialogTitle id='form-dialog-title'>
        Query History
        <Grid container className={classes.query_header}>
          <Grid xs={6} item>
            Comment
          </Grid>
          <Grid xs={4} item>
            Created
          </Grid>
        </Grid>
      </DialogTitle>
      <DialogContent>
        {savedQueries &&
          savedQueries.sort((a,b) => b.created_at.localeCompare(a.created_at)).map((query, i) => (
            <React.Fragment key={i}>
              <Grid container alignItems='center'>
                <Grid xs={6} item>
                  <Button variant='text' onClick={() => setIndex(selectedIndex == i ? -1 : i)}>
                    {query.comment || 'No comment'}
                  </Button>
                </Grid>
                <Grid xs={4} item>
                  <Typography>{formatTime(query.created_at)}</Typography>
                </Grid>
              </Grid>
              <Collapse in={selectedIndex == i} timeout='auto' unmountOnExit>
                <Grid container className={classes.diff}>
                  {selectedIndex == i && query.unpackedQuery && (
                    <Editor value={query.unpackedQuery}/>
                  )}
                </Grid>
              </Collapse>
            </React.Fragment>
          ))}
      </DialogContent>
      <DialogActions className={classes.actions}>
        {savedQueries && selectedIndex != -1 && (<>
          <Button onClick={() => removeQuery(savedQueries[selectedIndex])}>
            Remove
          </Button>
          <Button color='secondary' onClick={() => loadQuery(savedQueries[selectedIndex])}>
            Load
          </Button>
        </>)}
      </DialogActions>
    </Dialog>
  );
};

export default QueryHistory;
