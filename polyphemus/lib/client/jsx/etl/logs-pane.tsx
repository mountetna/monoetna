import React, {useEffect, useState, useMemo} from 'react';
import EtlPane, {EtlPaneHeader} from './etl-pane';
import {makeStyles} from '@material-ui/core/styles';
import {json_get} from 'etna-js/utils/fetch';

import {basicSetup, EditorView} from 'codemirror';
import {defaultHighlightStyle, syntaxHighlighting} from '@codemirror/language';
import {EditorState} from '@codemirror/state';
import {gutter, lineNumbers} from '@codemirror/view';
import CodeMirror from 'rodemirror';
import {lintGutter} from '@codemirror/lint';

import ErrorBoundary from 'etna-js/components/error_boundary';

const useStyles = makeStyles((theme) => ({
  editor: {
    border: '1px solid #ccc',
    height: 'calc(100vh - 375px)',
    marginBottom: '5px',
    resize: 'vertical',
    overflow: 'hidden',
    overflowY: 'auto'
  }
}));

const LogsPane = ({
  selected,
  config_id,
  name,
  project_name
}: {
  selected: string | null;
  config_id: number;
  name: string;
  project_name: string;
}) => {
  const classes = useStyles();
  const [output, setOutput] = useState('');

  useEffect(() => {
    json_get(`/api/etl/${project_name}/output/${config_id}`).then(({output}) =>
      setOutput(output)
    );
  }, []);

  const extensions = useMemo(
    () => [
      basicSetup,
      syntaxHighlighting(defaultHighlightStyle, {fallback: true}),
      EditorView.editable.of(false),
      EditorState.readOnly.of(true),
      EditorView.lineWrapping,
      EditorState.tabSize.of(2),
      gutter({class: 'CodeMirror-lint-markers'}),
      lineNumbers(),
      lintGutter()
    ],
    []
  );

  return (
    <EtlPane mode='logs' selected={selected}>
      <EtlPaneHeader title='Logs' />
      <div className={classes.editor}>
        <ErrorBoundary>
          <CodeMirror extensions={extensions} value={output} />
        </ErrorBoundary>
      </div>
      {output && (
        <a
          href={URL.createObjectURL(new Blob([output], {type: 'text/plain'}))}
          download={`log-${project_name}-${config_id}-${name}.txt`}
        >
          Download log
        </a>
      )}
    </EtlPane>
  );
};

export default LogsPane;
