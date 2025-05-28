import React, {useState, useEffect, useMemo} from 'react';
import {basicSetup, EditorView} from 'codemirror';
import {javascript} from '@codemirror/lang-javascript';
import {defaultHighlightStyle, syntaxHighlighting} from '@codemirror/language';
import {EditorState} from '@codemirror/state';
import {gutter, lineNumbers} from '@codemirror/view';
import CodeMirror from 'rodemirror';
import {lintGutter} from '@codemirror/lint';
import ErrorBoundary from './error_boundary';

const Editor = ({value}:{
  value: string
}) => {

  const extensions = useMemo(
    () => [
      basicSetup,
      syntaxHighlighting(defaultHighlightStyle, {fallback: true}),
      javascript(),
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

  return <ErrorBoundary>
    <CodeMirror
      extensions={extensions}
      value={value}
    />
  </ErrorBoundary>
}

export default Editor;
