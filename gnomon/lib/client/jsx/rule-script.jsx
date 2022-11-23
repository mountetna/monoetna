// Framework libraries.
import React, {useMemo} from 'react';

import {json, jsonParseLinter} from '@codemirror/lang-json';
import {defaultHighlightStyle, syntaxHighlighting} from '@codemirror/language';
import {basicSetup, EditorView} from 'codemirror';
import {EditorState} from '@codemirror/state';
import {gutter, lineNumbers} from '@codemirror/view';
import CodeMirror from 'rodemirror';
import {linter, lintGutter} from '@codemirror/lint';

import ErrorBoundary from 'etna-js/components/error_boundary';

function RuleScript({update,script}) {
  const extensions = [
    basicSetup,
    syntaxHighlighting(defaultHighlightStyle, {fallback: true}),
    json(),
    EditorView.editable.of(true),
    EditorState.tabSize.of(2),
    gutter({class: 'CodeMirror-lint-markers'}),
    EditorState.readOnly.of(false),
    EditorView.lineWrapping,
    linter(jsonParseLinter()),
    lineNumbers(),
    lintGutter(),
    EditorView.updateListener.of(function (e) {
      if (e.docChanged) {
        update(e.state.doc.toString());
      }
    })
  ];

  return (
    <div>
      <ErrorBoundary>
        <CodeMirror extensions={extensions} value={script} />
      </ErrorBoundary>
    </div>
  );
}
export default RuleScript;
