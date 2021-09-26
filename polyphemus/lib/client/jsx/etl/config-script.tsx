// Framework libraries.
import React, { useState, useEffect } from 'react';
import { Controlled } from 'react-codemirror2';

import CodeMirror, { Editor } from 'codemirror';
import {makeStyles} from '@material-ui/core/styles';

import 'codemirror/mode/javascript/javascript';
import 'codemirror/addon/lint/lint';

import JsonMap from 'json-source-map';

import Ajv from 'ajv';

const ajv = new Ajv({
  allErrors: true,
  verbose: true
});

const errorMessage = (error : any) => {
  let addendum;
  switch(error.keyword) {
    case 'additionalProperties':
      addendum = error.params.additionalProperty;
      break;
    case 'enum':
      addendum = error.params.allowedValues.join(', ');
      break;
  };
  return [
    error.instancePath,
    error.message,
    addendum
  ].filter(_=>_).join(' ');
}

export const validator = (schema:any, editor:Editor): Function => {
  const validate = ajv.compile(schema);

  return (text:string) => {
    // first see if the json can parse
    let json:any;

    // ignore blank texts
    if (!text) return [];

    try {
      json = JsonMap.parse(text);
    } catch (error) {
      // parse the error message and report the line numbers
      let { message } = error;

      let match = message.match(/^(?<message>[\s\S]*) at position (?<pos>\d+)/m);
      if (match) {
        return [
          {
            from: editor.posFromIndex(parseInt(match.groups.pos)),
            to: editor.posFromIndex(parseInt(match.groups.pos)),
            severity: 'error',
            message: match.groups.message
          }
        ]
      }
      return [];
    }

    let valid = validate(json.data);

    if (valid || !validate.errors) return [];

    return validate.errors!.map(
      error => {
        let pointer = json.pointers[error.instancePath];

        if (!pointer) return null;

        return {
          from: CodeMirror.Pos(pointer.value.line, pointer.value.column),
          to: CodeMirror.Pos(pointer.valueEnd.line, pointer.valueEnd.column),
          message: errorMessage(error),
          severity: 'error'
        }
      }
    ).filter(_=>_);
  }
}

const useStyles = makeStyles( theme => ({
  editor: {
    border: '1px solid #ccc',
    height: '200px'
  }
}));

const ConfigScript = ({ schema, script, onChange } : { schema: any, script: string, onChange: Function } ) => {
  let [ editor, setEditor ] = useState<Editor| null>(null);

  const classes = useStyles();

  useEffect(
    () => {
      if (editor) CodeMirror.registerHelper("lint", "json", validator(schema, editor))
    }, [ editor ]
  );

  return (
    <div className={classes.editor}>
      <Controlled
        options = {{
          readOnly: false,
          lineNumbers: true,
          lineWrapping: true,
          mode: 'application/json',
          gutters: ['CodeMirror-lint-markers'],
          lint: true,
          tabSize: 2
        }}
        editorDidMount={ setEditor }
        value={script}
        onBeforeChange={(editor, data, value) => { onChange(value); }}
      />
    </div>
  );
};

export default ConfigScript;
