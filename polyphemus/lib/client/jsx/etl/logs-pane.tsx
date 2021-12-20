import React, { useEffect, useState, useCallback } from 'react';
import Button from '@material-ui/core/Button';
import EtlPane, {EtlPaneHeader} from './etl-pane';
import {makeStyles} from '@material-ui/core/styles';
import { Controlled } from 'react-codemirror2';
import Typography from '@material-ui/core/Typography';
import { json_get } from 'etna-js/utils/fetch';

const useStyles = makeStyles( theme => ({
  editor: {
    border: '1px solid #ccc',
    height: '200px',
    resize: 'vertical',
    overflow: 'hidden'
  }
}));

const LogsPane = ({selected, name, project_name}:{selected:string|null,name:string,project_name:string}) => {
  const classes = useStyles();
  const [ output, setOutput ] = useState('');

  useEffect( () => {
    json_get(
      `/api/etl/${project_name}/output/${name}`
    ).then(
      ({output}) => setOutput(output)
    )
  }, []);

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
