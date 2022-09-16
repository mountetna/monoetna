// Framework libraries.
import React, {useMemo} from 'react';

import {json, jsonParseLinter} from '@codemirror/lang-json';
import {defaultHighlightStyle, syntaxHighlighting} from '@codemirror/language';
import {basicSetup, EditorView} from 'codemirror';
import {EditorState} from '@codemirror/state';
import {gutter, lineNumbers} from '@codemirror/view';
import CodeMirror from 'rodemirror';
import {linter, lintGutter} from '@codemirror/lint';

import ErrorBoundary from '../query/error_boundary';

function ViewScript(props) {
  let {is_editing, onChange, script} = props;
  let className = is_editing
    ? 'manifest-script editing view-editor'
    : 'manifest-script view-editor';

  const extensions = useMemo(() => {
    let config = [
      basicSetup,
      syntaxHighlighting(defaultHighlightStyle, {fallback: true}),
      json(),
      EditorView.editable.of(is_editing),
      EditorState.tabSize.of(2),
      gutter({class: 'CodeMirror-lint-markers'}),
      EditorState.readOnly.of(is_editing ? false : 'no-cursor'),
      EditorView.lineWrapping,
      linter(jsonParseLinter()),
      EditorView.updateListener.of(function (e) {
        if (e.docChanged) {
          onChange(e.state.doc.toString());
        }
      })
    ];

    if (is_editing) {
      config.push(lineNumbers());
      config.push(lintGutter());
    }

    return config;
  }, [is_editing]);

  return (
    <div className={className}>
      <ErrorBoundary>
        <CodeMirror extensions={extensions} value={script} />
      </ErrorBoundary>
    </div>
  );
}
export default ViewScript;
