import React, {useState, useEffect, useCallback, useContext} from 'react';

import Grid from '@material-ui/core/Grid';
import Link from '@material-ui/core/Link';
import TextField from '@material-ui/core/TextField';
import Button from '@material-ui/core/Button';
import IconButton from '@material-ui/core/IconButton';
import PlusOneIcon from '@material-ui/icons/PlusOne';
import FormControl from '@material-ui/core/FormControl';
import InputLabel from '@material-ui/core/InputLabel';
import FormHelperText from '@material-ui/core/FormHelperText';
import Select from '@material-ui/core/Select';
import MenuItem from '@material-ui/core/MenuItem';
import Typography from '@material-ui/core/Typography';
import ProjectHeader from 'etna-js/components/project-header';
import {makeStyles} from '@material-ui/core/styles';
import { json_get, json_post } from 'etna-js/utils/fetch';
import { magmaPath } from 'etna-js/api/magma_api';

import Letter from './letter';
import Bracket from './bracket';
import Corner from './corner';
import Line from './line';

import {isAdmin} from 'etna-js/utils/janus';

const useStyles = makeStyles((theme) => ({
  header: {
    borderBottom: '1px solid #eee'
  },
  tokens: {
    width: '100%',
    position: 'absolute'
  },
  composer: {
    marginLeft: 20,
    height: 'calc(100vh - 61px - 48px)',
    position: 'relative',
    overflowY: 'clip'
  },
  resolved: {
  },
  unresolved: {
    background: '#aaa',
    color: '#fff',
    borderBottom: '1px solid black'
  },
  counter: {
    background: '#fff',
    color: '#aaa',
    borderBottom: '1px solid black'
  },
  token: {
    fontSize: '30px'
  }
}));

const Letters = ({seq, className }) => {
  return <Grid container style={{ width: 'auto' }} className={className}>
    {
      seq.split('').map(
        (l,i) => <Letter key={i} letter={l}/>
      )
    }
  </Grid>
}

const firstValue = token => Object.values(token.values)[0];
const firstKey = token => Object.keys(token.values)[0];

const tokenEditorStyles = makeStyles((theme) => ({
  token_editor: {
    '&:hover $pointer, &:focus $pointer': {
      opacity: 1
    }
  },
  pointer: {
    position: 'absolute',
    opacity: 0.1,
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

const ResolvedEditor = ({token}) => (
  <TextField
    label={token.label}
    InputLabelProps={{ shrink: true }}
    InputProps={{ readOnly: true }}
    value={ firstValue(token) }
  />
);

const CounterEditor = ({token, tokens, value, update, pos, seq, project_name}) => (
  <React.Fragment>
    <TextField
      onChange={ e => update(pos, e.target.value) }
      value={value}
      InputLabelProps={{ shrink: true }}
      size='small'
      inputProps={{ pattern: "[0-9]*" }}
      label={token.label}/>
    {
      !value && token.filled && <IconButton
        onClick={ () => {
          const pre = tokens.slice(0,pos).map( t => t.seq ).join('');
          json_post(magmaPath(`gnomon/${project_name}/increment/${token.label.replace('_counter','')}/${pre}`)).then(
            value => update(pos, value)
          )
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

const TokenEditor = params => {
  const {token, seq, height, value, pos, update} = params;
  const classes = tokenEditorStyles()

  if (token.type == 'hidden') return null;

  const EditorComponent = eval( `${ token.type.charAt(0).toUpperCase() + token.type.slice(1) }Editor`);

  const voff = 40;
  const lh = 70;

  const w = 75;

  return <Grid className={classes.token_editor} style={{position: 'absolute', left: 0, top: 0}}>
    <Bracket
      bottom={voff}
      left={token.from * 40+2}
      width={(token.to - token.from + 1) * 40 - 4}
    />
    <Grid container alignItems='center' className={classes.editor}
      style={{
        left: w + 25 + seq.length * 40,
        bottom: (height / 2 - token.height) * lh,
        background: (token.type == 'resolved' || value) ? 'none' : '#eee'
      }}>
      <EditorComponent classes={classes} {...params}/>
    </Grid>
    <Corner
      className={ classes.pointer }
      left={ (token.to + token.from + 1) * 20 }
      bottom={ voff + 12.5 }
      right={ seq.length * 40 }
      top={ voff + 12 + 15 * (height - token.height + 1) } />
    <Line 
      className={ classes.pointer }
      left={ seq.length * 40 }
      right={ seq.length * 40 + w }
      top={ 7 + (height / 2 - token.height) * lh }
      bottom={ voff + 12 + 15 * (height - token.height + 1) }
    />
    <Line
      className={ classes.pointer }
      left={ seq.length * 40 + w }
      right={ seq.length * 40 + w + 20 }
      bottom={ 7 + (height / 2 - token.height) * lh }
      top={ 7 + (height / 2 - token.height) * lh }
    />
  </Grid>
}

const Token = ({token, value}) => {
  const classes = useStyles();

  return <Grid style={{position: 'relative'}}>
    <Letters className={classes[(value ? 'resolved' : token.type)]} seq={ token.seq }/>
  </Grid>
}

const ComposeIdentifier = ({project_name, rule_name}) => {
  const classes = useStyles();

  // 'MVIR1-HS169-D0PL1-CTK1';

  // a string of tokens we must satisfy
  const [ tokens, setTokens ] = useState([]);
  const [ values, setValues ] = useState([]);
  
  useEffect( () => {
    json_get(magmaPath(`gnomon/${project_name}/rule/${rule_name}`)).then(
      ({rule}) => {
        setTokens(rule);
        setValues( rule.map(t => '') );
      }
    )
  }, [] );

  const setValue = useCallback( (i, val) => {
    let newValues = [ ... values ];
    newValues[i] = val;
    console.log({newValues});
    setValues( newValues );
  }, [ values ]);

  // compute some ordering, etc. information about tokens for display and write it to the token
  const updateToken = useCallback(([ pos, height, filled ], token, i) => {
    let seq, type, new_filled;

    // there is only one option
    if (token.values && Object.keys(token.values).length == 1) {
      seq = firstKey(token);
      type = (firstValue(token)[0] == '#') ? 'hidden' : 'resolved'; 
      new_filled = filled;
    } else {
      seq = values[i] || token.name;
      type = (token.name == 'n') ? 'counter' : 'unresolved';
      new_filled = filled && !!values[i]
    }

    if (type != 'hidden') height = height + 1;

    Object.assign(token, { seq, type, height, from: pos, to: pos + seq.length - 1, filled });

    return [ pos + seq.length, height, new_filled ]
  }, [values]);

  const [ _, height ] = tokens.reduce( updateToken, [ 0, 0, true ] );

  const seq = tokens.map( t => t.seq ).join('');

  return <Grid>
    <ProjectHeader project_name={ project_name } className={classes.header}/>
    <Grid container alignItems='center' className={classes.composer} style={{ width: 40 * (seq.length+1) }}>
      <Grid container className={classes.tokens}>
      {
        tokens.map(
          (token, i) => <TokenEditor key={i} token={token} seq={seq} tokens={tokens} height={height} update={setValue} value={ values[i] } pos={i} project_name={project_name}/>
        )
      }
      </Grid>
      {
        tokens.map(
          (token, i) => <Token key={i} token={token} value={ values[i] }/>
        )
      }
    </Grid>
  </Grid>
}

export default ComposeIdentifier;
