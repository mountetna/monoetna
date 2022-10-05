// Framework libraries.
import React, {useMemo} from 'react';
import {
  defaultHighlightStyle,
  syntaxHighlighting,
  StreamLanguage
} from '@codemirror/language';
import {basicSetup, EditorView} from 'codemirror';
import {EditorState} from '@codemirror/state';
import {gutter, lineNumbers} from '@codemirror/view';
import CodeMirror from 'rodemirror';
import {simpleMode} from '@codemirror/legacy-modes/mode/simple-mode';

import ErrorBoundary from 'etna-js/components/error_boundary';

const timur_lang = [
  {regex: /\#.*$/m, token: 'comment'},
  {regex: /@[\w]+(?=\()/i, token: 'macro'},
  {regex: /@[\w]+/, token: 'variable'},
  {regex: /("|')(?:\\(?:\r\n|[\s\S])|(?!\1)[^\\\r\n])*\1/, token: 'string'},
  {regex: /[\w]+(?=\s*:)/i, token: 'label'},
  {regex: /\{/, token: 'template', next: 'template'},
  {regex: /[\w]+(?=\()/i, token: 'function'},
  {regex: /[\[\]]/, token: 'vector'},
  {regex: /\$[\w]+/, token: 'column'}
];

const template = [
  {regex: /\}/, token: 'template', next: 'start'},
  {regex: /\%[0-9]+/, token: 'template-var'},
  {regex: /[^}]/, token: 'template-text'}
];

const ManifestScript = ({is_editing, onChange, script}) => {
  let className = is_editing ? 'manifest-script editing' : 'manifest-script';

  const extensions = useMemo(() => {
    let config = [
      basicSetup,
      syntaxHighlighting(defaultHighlightStyle, {fallback: true}),
      EditorView.editable.of(is_editing),
      EditorState.tabSize.of(2),
      gutter({class: 'CodeMirror-lint-markers'}),
      EditorState.readOnly.of(is_editing ? false : 'no-cursor'),
      EditorView.lineWrapping,
      EditorView.updateListener.of(function (e) {
        if (e.docChanged) {
          onChange(e.state.doc.toString());
        }
      }),
      StreamLanguage.define(
        simpleMode({
          start: timur_lang,
          template,
          meta: {
            lineComment: '#'
          }
        })
      )
    ];

    if (is_editing) {
      config.push(lineNumbers());
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
};

export default ManifestScript;
