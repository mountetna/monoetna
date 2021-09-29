import React, { useState, useCallback } from 'react';
import Button from '@material-ui/core/Button';
import EtlPane, {EtlPaneHeader} from './etl-pane';
import {makeStyles} from '@material-ui/core/styles';
import { Controlled } from 'react-codemirror2';
import Typography from '@material-ui/core/Typography';

const useStyles = makeStyles( theme => ({
  editor: {
    border: '1px solid #ccc',
    height: '200px'
  }
}));

const LogsPane = ({selected, output}:{selected:string|null,output:string}) => {
  const classes = useStyles();

  return <EtlPane mode='logs' selected={selected}>
    <EtlPaneHeader title='Logs'/>
    <div className={classes.editor}>
      <Controlled
        options = {{
          readOnly: true,
          lineNumbers: true,
          lineWrapping: true,
          mode: 'text',
          gutters: ['CodeMirror-lint-markers'],
          lint: false,
          tabSize: 2
        }}
        value={output}
        onBeforeChange={(editor, data, value) => { }}
      />
    </div>
  </EtlPane>
}

export default LogsPane;
