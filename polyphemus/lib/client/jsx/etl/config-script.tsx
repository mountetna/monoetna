// Framework libraries.
import React, {useMemo} from 'react';
import {makeStyles} from '@material-ui/core/styles';

import {basicSetup, EditorView} from 'codemirror';
import {json} from '@codemirror/lang-json';
import {defaultHighlightStyle, syntaxHighlighting} from '@codemirror/language';
import {EditorState, Text} from '@codemirror/state';
import {gutter, lineNumbers} from '@codemirror/view';
import CodeMirror from 'rodemirror';
import {lintGutter, linter, Diagnostic, LintSource} from '@codemirror/lint';

import ErrorBoundary from 'etna-js/components/error_boundary';

import JsonMap from 'json-source-map';

import Ajv from 'ajv';

const ajv = new Ajv({
  allErrors: true,
  verbose: true
});

const errorMessage = (error: any) => {
  let addendum;
  switch (error.keyword) {
    case 'additionalProperties':
      addendum = error.params.additionalProperty;
      break;
    case 'enum':
      addendum = error.params.allowedValues.join(', ');
      break;
  }
  return [error.instancePath, error.message, addendum]
    .filter((_) => _)
    .join(' ');
};

function getErrorPosition(error: SyntaxError, doc: Text): number {
  let m;
  if ((m = error.message.match(/at position (\d+)/)))
    {return Math.min(+m[1], doc.length);}
  if ((m = error.message.match(/at line (\d+) column (\d+)/)))
    {return Math.min(doc.line(+m[1]).from + +m[2] - 1, doc.length);}
  return 0;
}

export const validator = (schema: any): LintSource => {
  const validate = ajv.compile(schema);

  return (view: EditorView): Diagnostic[] => {
    // first see if the json can parse
    let json: any;
    const text = view.state.doc.toString();

    // ignore blank texts
    if (!text) return [];

    try {
      json = JsonMap.parse(text);
    } catch (error) {
      // parse the error message and report the line numbers
      const pos = getErrorPosition(error, view.state.doc);

      return [
        {
          from: pos,
          to: pos,
          severity: 'error',
          message: error.message
        }
      ];
    }

    let valid = validate(json.data);

    if (valid || !validate.errors) return [];
    console.log(validate.errors);
    return validate.errors.reduce((acc: Diagnostic[], error) => {
      let pointer = json.pointers[error.instancePath];

      if (pointer) {
        console.log('pointer', pointer, error);
        acc.push({
          from: pointer.value.pos,
          to: pointer.valueEnd.pos,
          message: errorMessage(error),
          severity: 'error'
        } as Diagnostic);
      }

      return acc;
    }, []);
  };
};

const useStyles = makeStyles((theme) => ({
  editor: {
    border: '1px solid #ccc',
    height: '200px',
    resize: 'vertical',
    overflow: 'hidden',
    overflowY: 'auto'
  }
}));

const ConfigScript = ({
  schema,
  script,
  onChange
}: {
  schema: any;
  script: string;
  onChange: Function;
}) => {
  const classes = useStyles();

  const extensions = useMemo(
    () => [
      basicSetup,
      syntaxHighlighting(defaultHighlightStyle, {fallback: true}),
      json(),
      EditorView.editable.of(true),
      EditorState.readOnly.of(false),
      EditorView.lineWrapping,
      EditorState.tabSize.of(2),
      gutter({class: 'CodeMirror-lint-markers'}),
      lineNumbers(),
      lintGutter(),
      linter(validator(schema)),
      EditorView.updateListener.of(function (e) {
        if (e.docChanged) {
          onChange(e.state.doc.toString());
        }
      })
    ],
    []
  );

  return (
    <div className={classes.editor}>
      <ErrorBoundary>
        <CodeMirror extensions={extensions} value={script} />
      </ErrorBoundary>
    </div>
  );
};

export default ConfigScript;
