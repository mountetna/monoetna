import React, {useState, useContext, useMemo} from 'react';
import Grid from '@mui/material/Grid';
import IconButton from '@mui/material/IconButton';
import ExpandMoreIcon from '@mui/icons-material/ExpandMore';
import ExpandLessIcon from '@mui/icons-material/ExpandLess';
import Typography from '@mui/material/Typography';
import { makeStyles } from '@mui/styles';

import {QueryResultsContext} from '../../contexts/query/query_results_context';

import ErrorBoundary from 'etna-js/components/error_boundary';
import {json} from '@codemirror/lang-json';
import {defaultHighlightStyle, syntaxHighlighting} from '@codemirror/language';
import {EditorView} from 'codemirror';
import {EditorState} from '@codemirror/state';
import CodeMirror from 'rodemirror';

const useStyles = makeStyles((theme) => ({
  query_string: {
    width: '100%',
    padding: '10px',
    position: 'relative',
    borderBottom: '1px solid #eee'
  },
  query_box: {
    border: '1px solid #ccc',
    borderRadius: 2,
    backgroundColor: 'rgba(128,128,128,0.01)'
  },
  expand: {
    position: 'absolute',
    right: '10px',
    top: '10px',
    zIndex: '2'
  }
}));

const QueryString = () => {
  const classes = useStyles();
  const {state: { queryString }} = useContext(QueryResultsContext);

  const [ fold, setFold ] = useState(true)

  const extensions = useMemo(
    () => [
      syntaxHighlighting(defaultHighlightStyle, {fallback: true}),
      json(),
      EditorView.editable.of(false),
      EditorState.readOnly.of(true),
      EditorView.lineWrapping
    ],
    []
  );

  if (!queryString || queryString == '""') return null;

  const prettyString = JSON.stringify(JSON.parse(queryString), null, 2);

  const codeMirrorText = fold ? queryString : prettyString;

  return (
    <Grid className={ classes.query_string } >
      <IconButton className={ classes.expand } onClick={ () => setFold(!fold) } >
        { fold ? <ExpandMoreIcon/> : <ExpandLessIcon/> }
      </IconButton>
      <Grid className={classes.query_box}>
        <ErrorBoundary>
          <CodeMirror extensions={extensions} value={codeMirrorText} />
        </ErrorBoundary>
      </Grid>
    </Grid>
  );
};

export default QueryString;
