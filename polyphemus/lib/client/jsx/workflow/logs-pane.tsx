import React, {useEffect, useState, useMemo, useCallback} from 'react';
import WorkflowPane, {WorkflowPaneHeader} from './workflow-pane';

import Collapse from '@material-ui/core/Collapse';
import Table from '@material-ui/core/Table';
import TableBody from '@material-ui/core/TableBody';
import TableCell from '@material-ui/core/TableCell';
import TableSortLabel from '@material-ui/core/TableSortLabel';
import TableRow from '@material-ui/core/TableRow';
import TableContainer from '@material-ui/core/TableContainer';
import TableHead from '@material-ui/core/TableHead';
import {makeStyles} from '@material-ui/core/styles';
import {json_get} from 'etna-js/utils/fetch';

import {basicSetup, EditorView} from 'codemirror';
import {defaultHighlightStyle, syntaxHighlighting} from '@codemirror/language';
import {EditorState} from '@codemirror/state';
import {gutter, lineNumbers} from '@codemirror/view';
import CodeMirror from 'rodemirror';
import {lintGutter} from '@codemirror/lint';
import {Run} from '../polyphemus';
import {formatTime} from './run-state';

import ErrorBoundary from 'etna-js/components/error_boundary';

const useStyles = makeStyles((theme) => ({
  run_list: {
    border: '1px solid #ccc',
    height: 'calc(100vh - 330px)'
  },
  run: {
    cursor: 'pointer',
    '&:hover': {
      background: '#eee'
    }
  },
  output_view: {
    padding: '0px',
    border: 'none'
  },
  editor: {
    borderBottom: '1px solid #ccc',
    height: '500px',
    marginBottom: '5px',
    resize: 'vertical',
    overflow: 'hidden',
    overflowY: 'auto'
  }
}));

const extensions = [
  basicSetup,
  syntaxHighlighting(defaultHighlightStyle, {fallback: true}),
  EditorView.editable.of(false),
  EditorState.readOnly.of(true),
  EditorView.lineWrapping,
  EditorState.tabSize.of(2),
  gutter({class: 'CodeMirror-lint-markers'}),
  lineNumbers(),
  lintGutter()
];

const RunRow = ({run, project_name} : { run: Run; project_name: string; }) => {
  const classes = useStyles();
  const { run_id, created_at, config_id, version_number, status, finished_at } = run;
  const [ showOutput, setShowOutput ] = useState(false);
  const [ output, setOutput ] = useState<string|null>(null);
  const toggleOutput = useCallback(() => {
    if (showOutput) {
      setShowOutput(false);
      return;
    }

    json_get(`/api/workflows/${project_name}/run/output/${run_id}`).then(({output}) => {
      setShowOutput(true);
      setOutput(output);
    });
  }, [showOutput]);
  return <>
    <TableRow className={classes.run} onClick={ toggleOutput }>
      <TableCell>{run_id}</TableCell>
      <TableCell>{version_number}</TableCell>
      <TableCell>{status}</TableCell>
      <TableCell>{formatTime(created_at)}</TableCell>
      <TableCell>{finished_at ? formatTime(finished_at) : 'never run'}</TableCell>
    </TableRow>
    <TableRow>
      <TableCell colSpan={5} className={classes.output_view}>
        <Collapse in={ showOutput }>
          <div className={classes.editor}>
            <ErrorBoundary>
              <CodeMirror extensions={extensions} value={output || ''} />
            </ErrorBoundary>
          </div>
          {output && (
            <a
              href={URL.createObjectURL(new Blob([output], {type: 'text/plain'}))}
              download={`log-${project_name}-${config_id}-${run_id}.txt`}
            >
              Download log
            </a>
          )}
        </Collapse>
      </TableCell>
    </TableRow>
  </>
}

const LogsPane = ({
  selected,
  config_id,
  project_name
}: {
  selected: string | null;
  config_id: number;
  project_name: string;
}) => {
  const classes = useStyles();
  const [output, setOutput] = useState('');
  const [runs, setRuns] = useState<Run[]>([]);

  useEffect(() => {
    json_get(`/api/workflows/${project_name}/runs/${config_id}`).then(({runs}) =>
      setRuns(runs)
    );
  }, []);


  return (
    <WorkflowPane mode='logs' selected={selected}>
      <WorkflowPaneHeader title='Runs' />
      <TableContainer className={classes.run_list}>
        <Table size="small" stickyHeader>
          <TableHead>
            <TableRow>
              <TableCell>Run id</TableCell>
              <TableCell>Config Version</TableCell>
              <TableCell>Status</TableCell>
              <TableCell>Created</TableCell>
              <TableCell>Finished</TableCell>
            </TableRow>
          </TableHead>
          <TableBody>
          {
            runs.sort( (a,b) => a.created_at.localeCompare(b.created_at)).map(
              run => <RunRow project_name={project_name} key={run.run_id} run={run}/>
            )
          }
          </TableBody>
        </Table>
      </TableContainer>
    </WorkflowPane>
  );
};

export default LogsPane;
