import React, { useState, useEffect, useCallback } from 'react';
import {makeStyles} from '@material-ui/core/styles';
import Grid from '@material-ui/core/Grid';
import Button from '@material-ui/core/Button';
import Dialog from '@material-ui/core/Dialog';
import DialogActions from '@material-ui/core/DialogActions';
import DialogContent from '@material-ui/core/DialogContent';
import DialogTitle from '@material-ui/core/DialogTitle';

import Collapse from '@material-ui/core/Collapse';
import { json_get } from 'etna-js/utils/fetch';
import { formatTime } from './run-state';
import { Controlled } from 'react-codemirror2';

import { createTwoFilesPatch } from 'diff';

import 'codemirror/mode/javascript/javascript';
import 'codemirror/mode/diff/diff';
import 'codemirror/addon/lint/lint';

const useStyles = makeStyles( theme => ({
  revision_header: {
    borderBottom: '1px solid #ccc',
    padding: '5px 0px'
  },
  diff: {
    height: '200px',
    border: '1px solid #ccc'
  },
  actions: {
    paddingRight: '25px'
  }
}));

const diff = ({revision, prev, diffType, config}:{
  revision:EtlRevision,
  prev:EtlRevision,
  diffType:string,
  config:any
}) => {
  if (diffType == 'text') return JSON.stringify(revision.config,null,2);

  if (diffType == 'curr')
    return createTwoFilesPatch(
      'revision',
      'curr',
      JSON.stringify(revision.config,null,2)+'\n',
      JSON.stringify(config,null,2)+'\n'
    );
  else
    return createTwoFilesPatch(
      'prev',
      'revision',
      JSON.stringify(prev ? prev.config : {},null,2)+'\n',
      JSON.stringify(revision.config,null,2)+'\n'
    );
}

const RevisionHistory = ({project_name, name, update, open, onClose, config}:{
  project_name: string,
  name: string,
  update: (config: any) => void,
  open: boolean,
  onClose: (event:{})=>void,
  config: any
}) => {
  const [ revisions, setRevisions ] = useState<null | EtlRevision[]>(null);
  const [ selectedRevision, setSelectedRevision ] = useState<number|null>(null);
  const [ diffType, setDiffType ] = useState('text');

  const setSelected = (revision:number, diffType:string) => {
    setSelectedRevision(revision);
    setDiffType(diffType);
  }

  const classes = useStyles();

  useEffect( () => {
    json_get(
      `/api/etl/${project_name}/revisions/${name}`
    ).then(
      revisions => setRevisions(revisions)
    )
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
              <Grid xs={4} item>{formatTime(revision.updated_at) }</Grid>
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
                  value={diff({revision, diffType, config, prev: revisions[i+1]})}
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
        revisions && selectedRevision != null && <Button onClick={() => update(revisions[selectedRevision].config)} >Load</Button>
      }
    </DialogActions>
  </Dialog>;
}

export default RevisionHistory;
