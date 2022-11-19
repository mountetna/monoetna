import React, {useState, useReducer, useEffect, useCallback, useContext, useMemo} from 'react';

import Grid from '@material-ui/core/Grid';
import Link from '@material-ui/core/Link';
import TextField from '@material-ui/core/TextField';
import Button from '@material-ui/core/Button';
import IconButton from '@material-ui/core/IconButton';
import ProjectHeader from 'etna-js/components/project-header';
import {makeStyles} from '@material-ui/core/styles';
import { json_get, json_post } from 'etna-js/utils/fetch';
import { magmaPath, getDocuments } from 'etna-js/api/magma_api';

import Letter from './letter';
import TokenEditor from './token-editor';

import Table from '@material-ui/core/Table';
import TableBody from '@material-ui/core/TableBody';
import TableCell from '@material-ui/core/TableCell';
import TableContainer from '@material-ui/core/TableContainer';
import TableHead from '@material-ui/core/TableHead';
import TableRow from '@material-ui/core/TableRow';
import Paper from '@material-ui/core/Paper';
import { TableWithTitle, IdTreeTable, MatchingNamesTable } from './match-tables';
import { Typography, Tooltip }from '@material-ui/core';

import {isAdmin} from 'etna-js/utils/janus';

require('../img/distort.svg');
require('../img/distort2.svg');
require('../img/distort3.svg');

const useStyles = makeStyles((theme) => ({
  header: {
    borderBottom: '1px solid #eee'
  },
  tables: {
    flex: '0 0 275px',
    width: '100%'
  },
  table_column: {
    flex: '0 0 50%',
    '&:first-of-type': {
      borderRight: '1px solid #ccc'
    },
    borderTop: '1px solid #ccc',
    height: 275,
    overflowY: 'scroll'
  },
  tokens: {
    position: 'absolute'
  },
  mainContent:{
    height: 'calc(100vh - 61px - 48px)',
    flexDirection: 'column',
    overflowY: 'auto',
    overflowX: 'hidden'
  },
  composer: {
    height: 'calc(100vh - 61px - 48px - 275px)',
    flex: '1 1 auto',
    width: '100%',
    position: 'relative',
    overflowY: 'scroll',
    paddingLeft: '10px',
    marginRight: '10px'
  },
  create: {
    position: 'absolute',
    top: 0,
    left: 20,
    width: 400
  },
  create_button: {
    backgroundImage: 'url(/images/distort.svg)',
    '&:hover': {
      backgroundImage: 'url(/images/distort2.svg)',
      boxShadow: '0 0 20px #d95'
    },
    '&:active': {
      backgroundImage: 'url(/images/distort3.svg)',
      boxShadow: '0 0 2px #d95'
    },
    backgroundSize: 'cover',
    border: '2px solid darkgoldenrod',
    boxShadow: '0 0 10px #d95',
    color: 'white',
    fontWeight: 'bold'
  },
  match_names: {
    width: '100%',
    height: '240px',
    overflowX: 'hidden',
    overflowY: 'scroll'
  },
  match_rules: {
    width: '100%',
    height: '240px',
    overflowY: 'scroll'
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

const Token = ({token, value}) => {
  const classes = useStyles();

  return <Grid style={{position: 'relative'}}>
    <Letters className={classes[(value ? 'resolved' : token.type)]} seq={ token.seq }/>
  </Grid>
}

function matchIds(ids, idRegex) {
  if (ids == null) return null
  const re = new RegExp('^'+idRegex+'$')
  return ids.filter(id => re.test(id.identifier))
}

const CreateButton = ({project_name, rule_name, identifier, onMouseOver, onMouseOut, update}) => {
  const classes = useStyles();
  return <Button
    onClick={ update }
    className={classes.create_button}
    onMouseOver={onMouseOver}
    onMouseOut={onMouseOut}
    size='large'
    color='primary'>Create</Button>
}

const firstValue = token => Object.values(token.values)[0];
const firstKey = token => Object.keys(token.values)[0];

const isSingleValueToken = token => token.values && Object.keys(token.values).length == 1;
const isHiddenToken = token => firstValue(token)[0] == '#';
const isCounterToken = token => token.name == 'n';

const regexForTokens = (tokens) => {
  if (tokens == null || tokens == []) return null;

  const regex = tokens.map( (token, i) => {
    if (['resolved', 'hidden'].includes(token.type)) return token.seq;
    else if (token.value != '') return token.value;
    else if (token.type == 'counter') return '[0-9]+';
    else return `(${ Object.keys(token.values).join('|') })`;
  }).join('');

  return regex;
}

const setTokens = (state, action) => {
  const [ _, tokens ] = action.tokens.reduce(
    // - 'filled' means all *previous* tokens are assigned, viz., either a
    //   single-value token or have a user-specified value
    // - 'assigned' means your current token has an assigned value
    // - 'type' is a fixed class for displaying the token and does not change
    // - 'seq' is the current displayed string (either the token placeholder or a value)
    // - 'height' is the cumulative count of visible tokens (mostly ignoring separators)
    ([filled, newTokens], token, i) => {
      let seq, type, assigned;
      let new_filled = filled;
      let height = i ? newTokens[i-1].height : 0;
      let from = i ? newTokens[i-1].to : 0;
      if (isSingleValueToken(token)) {
        seq = firstKey(token);
        type = isHiddenToken(token) ? 'hidden' : 'resolved';
        assigned = true;
      } else {
        seq = token.name;
        type = isCounterToken(token) ? 'counter' : 'unresolved';
        new_filled = false;
        assigned = false;
      }

      if (type != 'hidden') height = height + 1;

      newTokens.push({ ...token, seq, filled, assigned, type, height, value: '', from, to: from + seq.length });

      return [ new_filled, newTokens ];
    }, [ true, [] ]
  );

  return {
    ...state,
    tokens,
    seq: tokens.map( t => t.seq ).join(''),
    height: tokens.length ? tokens[tokens.length - 1].height : 0,
    regex: regexForTokens(tokens)
  }
}

const setValue = (state, action) => {
  const token = state.tokens[action.pos];
  const seq = action.value || token.name;
  const prevToken = action.pos ? state.tokens[action.pos - 1] : null;
  const newToken = {
    ...token,
    seq,
    filled: !prevToken || (prevToken.filled && prevToken.assigned),
    to: token.from + seq.length,
    value: action.value,
    assigned: !!action.value
  };

  const tokens = [
    ...state.tokens.slice(0,action.pos),
    newToken,
    ...state.tokens.slice(action.pos+1).reduce(
      (afterTokens, token, i) => {
        const from = i ? afterTokens[i-1].to : newToken.to;
        afterTokens.push(
          {
            ...token,
            filled: i ? (afterTokens[i-1].filled && afterTokens[i-1].assigned) : (newToken.filled && newToken.assigned),
            from,
            to: from + token.seq.length
          }
        )
        return afterTokens;
      }, []
    )
  ];

  return {
    ...state,
    tokens,
    seq: tokens.map(t => t.seq).join(''),
    regex: regexForTokens(tokens)
  }
}
export const reducer = (state, action) => {
  switch(action.type) {
    case 'SET_TOKENS':
      return setTokens(state, action);
    case 'SET_VALUE':
      return setValue(state, action);
    default:
      return state;
  }
}

const ComposeIdentifier = ({project_name, rule_name}) => {
  // 'MVIR1-HS169-D0PL1-CTK1';
  const classes=useStyles();
  // a string of tokens we must satisfy
  const [ state, dispatch ] = useReducer(reducer, {
    tokens: [],
    seq: ''
  });

  const [ names, setNames ] = useState(null);
  const [ decomposition, setDecomposition ] = useState(null);
  const [ highlight, setHighlight ] = useState(false);
  
  const { seq, height, tokens, regex } = state;
  
  useEffect( () => {
    json_get(magmaPath(`gnomon/${project_name}/list/${rule_name}`)).then(
      id_list => setNames(id_list)
    )
    json_get(magmaPath(`gnomon/${project_name}/rule/${rule_name}`)).then(
      ({rule}) => dispatch({ type: 'SET_TOKENS', tokens: rule })
    )
  }, [] );

  useEffect( () => {
    if (tokens.every(({assigned}) => assigned)) {
      json_get(magmaPath(`gnomon/${project_name}/decompose/${regex}`)).then(
        decomposition => setDecomposition(decomposition)
      );
    } else {
      setDecomposition(null);
    }
  }, [ tokens ]);

  const lh = 70;
  const dh = Math.max(height * lh + 10, 150);
  const dw = Math.max(40 * (seq.length+1), 400);

  return <Grid>
    <ProjectHeader project_name={ project_name } className={classes.header}/>
    <Grid container direction='row' className={classes.mainContent}>
      <Grid item container className={classes.composer}>
        <Grid item container alignItems='center' style={{ height: dh }}>
          <Grid container className={classes.tokens} style={{ width: dw }}>
          {
            tokens.map(
              (token, i) => <TokenEditor
                key={i}
                lh={lh}
                token={token}
                seq={seq}
                tokens={tokens}
                height={height}
                update={ (pos, value) => dispatch({type: 'SET_VALUE', pos, value })}
                value={ token.value }
                pos={i}
                project_name={project_name}
              />
            )
          }
          </Grid>
          {
            tokens.map(
              (token, i) => <Token key={i} token={token} value={ token.value }/>
            )
          }
          <Grid item container direction='column' alignItems='center' className={classes.create}
            style={{ top: dh/2 + 90, left: dw / 2 - 162 }}>
            {
              decomposition && (!decomposition.rules[rule_name].name_created_at ? <>
                <CreateButton
                  project_name={project_name}
                  rule_name={rule_name}
                  identifier={seq}
                  update={
                    () => json_post(magmaPath(`gnomon/${project_name}/generate/${rule_name}/${seq}`))
                      .then( decomposition => setDecomposition(decomposition))
                      .catch( (e) => console.log(e) )
                  }
                  onMouseOver={ () => setHighlight(true) }
                  onMouseOut={ () => setHighlight(false) }
                  />
              </> : null)
            }
          </Grid>
        </Grid>
      </Grid>
      <Grid item container className={classes.tables}>
        <IdTreeTable
          title='Matching Rules'
          boxClassName={classes.table_column}
          tableClassName={classes.match_rules}
          decomposition={decomposition}
          project_name={project_name}
          markNotCreated={true}
          highlight={ highlight }/>
        <MatchingNamesTable
          tableClassName={classes.match_names}
          names={matchIds(names, regex)}
          rule_name={rule_name}
          decomposition={decomposition}
          title='Matching Names' boxClassName={classes.table_column}/>
      </Grid>
    </Grid>
  </Grid>
}

export default ComposeIdentifier;
