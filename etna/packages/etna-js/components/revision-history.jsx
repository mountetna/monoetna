import React, { useState, useEffect, useCallback } from 'react';
import {makeStyles} from '@material-ui/core/styles';
import Grid from '@material-ui/core/Grid';
import Typography from '@material-ui/core/Typography';
import Button from '@material-ui/core/Button';
import Dialog from '@material-ui/core/Dialog';
import DialogActions from '@material-ui/core/DialogActions';
import DialogContent from '@material-ui/core/DialogContent';
import DialogTitle from '@material-ui/core/DialogTitle';

import Collapse from '@material-ui/core/Collapse';
import { Controlled } from 'react-codemirror2';

import { createTwoFilesPatch } from 'diff';

import 'codemirror/mode/javascript/javascript';
import 'codemirror/mode/diff/diff';
import 'codemirror/addon/lint/lint';

const formatTime = time => Intl.DateTimeFormat('en-US', {
   year: 'numeric', month: 'short', day: 'numeric',
   hour: 'numeric', minute: 'numeric'
}).format(new Date(time))

const useStyles = makeStyles( theme => ({
  revision_header: {
    borderBottom: '1px solid #ccc',
    padding: '5px 0px'
  },
  diff: {
    height: '200px',
    border: '1px solid #ccc',
    '& .react-codemirror2,.CodeMirror': {
      height: '100%',
      width: '100%'
    }
  },
  actions: {
    paddingRight: '25px'
  }
}));

const toJson = doc => JSON.stringify(doc,null,2);

const diff = ({revision, current, prev, revisionDoc, diffType}) => {
  if (diffType == 'text') return revisionDoc(revision);

  if (diffType == 'curr')
    return createTwoFilesPatch(
      'revision',
      'curr',
      revisionDoc(revision),
      revisionDoc(current)
    );
  else
    return createTwoFilesPatch(
      'prev',
      'revision',
      prev ? revisionDoc(prev) : revisionDoc({}),
      revisionDoc(revision)
    );
}

const RevisionHistory = ({getRevisions, revisionDoc, update, open, onClose}) => {
  const [ revisions, setRevisions ] = useState(null);
  const [ selectedRevision, setSelectedRevision ] = useState(null);
  const [ diffType, setDiffType ] = useState('text');

  const setSelected = (revision, newDiffType) => {
    if (selectedRevision == revision && diffType == newDiffType) {
      setSelectedRevision(null);
      return;
    }
    setSelectedRevision(revision);
    setDiffType(newDiffType);
  }

  const classes = useStyles();

  useEffect( () => {
    getRevisions().then(revisions => setRevisions(revisions))
  }, []);

  return <Dialog maxWidth="md" fullWidth scroll='paper' open={open} onClose={ onClose } aria-labelledby="form-dialog-title">
    <DialogTitle id="form-dialog-title">Revision History
      <Grid container className={classes.revision_header}>
        <Grid xs={6} item>Comment</Grid>
        <Grid xs={2} item>Diff</Grid>
        <Grid xs={4} item>Updated</Grid>
      </Grid>
    </DialogTitle>
    <DialogContent>
      {
        revisions && revisions.map( (revision,i) =>
          <React.Fragment key={i}>
            <Grid container alignItems='center'>
              <Grid xs={6} item >
                <Button variant='text' onClick={ () => setSelected(i, 'text') }>
                  { revision.comment || 'No comment'}
                </Button>
              </Grid>
              <Grid xs={2}  container item >
                { i != 0 ? <Button onClick={ () => setSelected(i,'curr') } variant='text' size='small'>curr</Button> : <div style={{width: '64px', display: 'inline'}}/> }
                { i < revisions.length-1 && <Button onClick={ () => setSelected(i, 'prev') } variant='text' size='small'>prev</Button> }
              </Grid>
              <Grid xs={4} item><Typography>{formatTime(revision.updated_at) }</Typography></Grid>
            </Grid>
            <Collapse in={selectedRevision == i} timeout='auto' unmountOnExit>
              <Grid container className={classes.diff}>
                {selectedRevision == i && <Controlled
                  options = {{
                    readOnly: true,
                    lineNumbers: true,
                    lineWrapping: true,
                    lint: true,
                    mode: diffType == 'text' ? 'javascript' : 'diff',
                    gutters: ['CodeMirror-lint-markers'],
                    tabSize: 2
                  }}
                  value={diff({revision, diffType, revisionDoc, current: revisions[0], prev: revisions[i+1]})}
                  onBeforeChange={(editor, data, value) => { }}
                /> }
              </Grid>
            </Collapse>
          </React.Fragment>
        )
      }
    </DialogContent>
    <DialogActions className={classes.actions}>
      {
        revisions && selectedRevision != null && <Button onClick={() => update(revisions[selectedRevision])} >Load</Button>
      }
    </DialogActions>
  </Dialog>;
}

export default RevisionHistory;
