import React, {useState, useEffect, useCallback, useContext} from 'react';

import Grid from '@material-ui/core/Grid';
import TextField from '@material-ui/core/TextField';
import IconButton from '@material-ui/core/IconButton';
import PlusOneIcon from '@material-ui/icons/PlusOne';
import FormControl from '@material-ui/core/FormControl';
import InputLabel from '@material-ui/core/InputLabel';
import Select from '@material-ui/core/Select';
import MenuItem from '@material-ui/core/MenuItem';
import Typography from '@material-ui/core/Typography';
import {makeStyles} from '@material-ui/core/styles';
import { json_get, json_post } from 'etna-js/utils/fetch';
import { magmaPath } from 'etna-js/api/magma_api';

import Bracket from './bracket';
import Corner from './corner';
import Line from './line';

const resolvedStyle = makeStyles( theme => ({
  resolved: {
    '& .MuiInputBase-root.Mui-disabled': {
      color: 'black'
    },
    '& .MuiInputLabel-root.Mui-disabled': {
      color: 'rgba(0,0,0,0.5)'
    }
  }
}) );

const firstValue = token => Object.values(token.values)[0];

const ResolvedEditor = ({token, description}) => {
  const classes = resolvedStyle();
  return <TextField
    className={classes.resolved}
    disabled
    label={token.label}
    InputLabelProps={{ shrink: true }}
    InputProps={{ readOnly: true }}
    value={ description || firstValue(token) }
  />;
};

const CounterEditor = ({token, tokens, value, update, pos, seq, project_name}) => (
  <React.Fragment>
    <TextField
      onChange={ e => update(pos, e.target.value) }
      value={value}
      InputLabelProps={{ shrink: true }}
      size='small'
      inputProps={{ pattern: '[0-9]*' }}
      label={token.label}/>
    {
      !value && token.filled && <IconButton
        onClick={ () => {
          const pre = tokens.slice(0,pos).map( t => t.seq ).join('');
          json_post(magmaPath(`gnomon/${project_name}/increment/${token.label.replace('_counter','')}/${pre}`)).then(
            value => update(pos, value)
          );
        }}
        size='small'
        title={ `Fill the next available value for ${token.label}`}>
        <PlusOneIcon/>
      </IconButton>
    }
  </React.Fragment>
);

const UnresolvedEditor = ({ token, value, pos, update, classes}) => (
  <FormControl style={{ width: 230 }}>
      <InputLabel shrink>{token.label}</InputLabel>
      <Select onChange={ e => update(pos, e.target.value) } value={value} displayEmpty size='small'>
        <MenuItem value=''><em className={classes.empty}>None</em></MenuItem>
        {
          Object.keys(token.values).map(
            name => <MenuItem key={name} value={name}>{token.values[name]}</MenuItem>
          )
        }
      </Select>
  </FormControl>
);

const tokenEditorStyles = makeStyles((theme) => ({
  token_editor: {
    '&:hover $pointer, &:focus $pointer': {
      opacity: 1
    }
  },
  pointer: {
    position: 'absolute',
    opacity: 0.2,
    borderTop: '1px solid black'
  },
  editor: {
    width: 290,
    boxShadow: '0 0 5px rgba(0,0,0,0.05)',
    position: 'absolute',
    border: '1px solid #888',
    padding: '5px',
    whiteSpace: 'nowrap',
  },
  empty: {
    color: '#888'
  }
}));

const TokenEditor = params => {
  const {token, voff=40, seq, lh=70, w=75, height, value, pos, update} = params;
  const classes = tokenEditorStyles();

  if (token.type == 'hidden') return null;

  const EditorComponent = eval( `${ token.type.charAt(0).toUpperCase() + token.type.slice(1) }Editor`);

  const ch = (height / 2 - token.height) * lh;

  return <Grid className={classes.token_editor} style={{position: 'absolute', left: 0, top: 0}}>
    <Bracket
      className={ classes.pointer }
      bottom={voff}
      left={token.from * 40+2}
      width={(token.to - token.from) * 40 - 4}
    />
    <Grid container alignItems='center' className={classes.editor}
      style={{
        left: w + 25 + seq.length * 40,
        bottom: ch,
        background: (token.type == 'resolved' || value) ? 'none' : '#eee'
      }}>
      <EditorComponent classes={classes} {...params}/>
    </Grid>
    <Corner
      className={ classes.pointer }
      left={ (token.to + token.from) * 20 }
      bottom={ voff + 11 }
      right={ seq.length * 40 }
      top={ voff + 12 + 15 * (height - token.height + 1) } />
    <Line 
      className={ classes.pointer }
      left={ seq.length * 40 }
      right={ seq.length * 40 + w }
      top={ 7 + ch }
      bottom={ voff + 12 + 15 * (height - token.height + 1) }
    />
    <Line
      className={ classes.pointer }
      left={ seq.length * 40 + w }
      right={ seq.length * 40 + w + 20 }
      bottom={ 7 + ch }
      top={ 7 + ch }
    />
  </Grid>;
};

export default TokenEditor;
