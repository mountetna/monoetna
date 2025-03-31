import React, {useState, useContext, useMemo} from 'react';
import Grid from '@material-ui/core/Grid';
import IconButton from '@material-ui/core/IconButton';
import ExpandMoreIcon from '@material-ui/icons/ExpandMore';
import ExpandLessIcon from '@material-ui/icons/ExpandLess';
import Typography from '@material-ui/core/Typography';
import { makeStyles } from '@material-ui/core/styles';

import {QueryResultsContext} from '../../contexts/query/query_results_context';

import ErrorBoundary from 'etna-js/components/error_boundary';
import {json} from '@codemirror/lang-json';
import {defaultHighlightStyle, syntaxHighlighting} from '@codemirror/language';
import {EditorView} from 'codemirror';
import {EditorState} from '@codemirror/state';
import CodeMirror from 'rodemirror';

const useStyles = makeStyles((theme) => ({
  query_string: {
    width: 'calc(100% - 20px)',
    padding: '10px',
    position: 'relative',
    borderBottom: '1px solid #eee'
  },
  query_box: {
    border: '1px solid #ccc',
    borderRadius: 2,
    paddingRight: '26px',
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

  const FoldIcon = fold ? ExpandMoreIcon : ExpandLessIcon;

  return (
    <Grid className={ classes.query_string } >
      <IconButton className={ classes.expand } onClick={ () => setFold(!fold) } size='small'>
        <FoldIcon fontSize='small'/>
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
