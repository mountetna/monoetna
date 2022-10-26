import React, {useState, useEffect, useCallback, useContext, useMemo} from 'react';

import Grid from '@material-ui/core/Grid';
import Link from '@material-ui/core/Link';
import TextField from '@material-ui/core/TextField';
import Button from '@material-ui/core/Button';
import IconButton from '@material-ui/core/IconButton';
import ProjectHeader from 'etna-js/components/project-header';
import {makeStyles} from '@material-ui/core/styles';
import { json_get, json_post } from 'etna-js/utils/fetch';
import { dateFormat } from 'etna-js/utils/format';
import { magmaPath, getDocuments } from 'etna-js/api/magma_api';

import Letter from './letter';
import TokenEditor, { firstKey, firstValue } from './token-editor';

import Table from '@material-ui/core/Table';
import TableBody from '@material-ui/core/TableBody';
import TableCell from '@material-ui/core/TableCell';
import TableContainer from '@material-ui/core/TableContainer';
import TableHead from '@material-ui/core/TableHead';
import TableRow from '@material-ui/core/TableRow';
import Paper from '@material-ui/core/Paper';

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
    height: 'calc(100vh - 61px - 48px - 250px)',
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
  ids_table: {
    height: 250,
    overflowX: 'hidden'
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
  const re = new RegExp("^"+idRegex+"$")
  return ids.filter((val) => re.test(val.identifier))
}

const IdsShower = ({ids}) => {
  const classes = useStyles();
  console.log('matches:', ids)
  console.log(ids == {})
  if (ids == null) return null
  if (ids.length === 0) return <div> No Matching Identifiers </div>
  return <TableContainer className={classes.ids_table} component={Paper}>
    <Table stickyHeader size="small">
      <TableHead>
        <TableRow>
          <TableCell>Name</TableCell>
          <TableCell align="right">Author</TableCell>
          <TableCell align="right">Named at</TableCell>
          <TableCell align="right">Recorded at</TableCell>
        </TableRow>
      </TableHead>
      <TableBody>
        {ids.map((id) => (
          <TableRow key={id.identifier}>
            <TableCell component="th" scope="row">
              {id.identifier}
            </TableCell>
            <TableCell align="right">{id.author}</TableCell>
            <TableCell align="right">{dateFormat(id.name_created_at)}</TableCell>
            <TableCell align="right">{dateFormat(id.record_created_at)}</TableCell>
          </TableRow>
        ))}
      </TableBody>
    </Table>
  </TableContainer>
}

const ComposeIdentifier = ({project_name, rule_name}) => {
  const classes = useStyles();

  // 'MVIR1-HS169-D0PL1-CTK1';

  // a string of tokens we must satisfy
  const [ tokens, setTokens ] = useState([]);
  const [ values, setValues ] = useState([]);

  const [ ids, setIds ] = useState(null)
  
  useEffect( () => {
    json_get(magmaPath(`gnomon/${project_name}/rule/${rule_name}`)).then(
      ({rule}) => {
        setTokens(rule);
        setValues( rule.map(t => '') );
      }
    )
  }, [] );

  useEffect( () => {
    json_get(magmaPath(`gnomon/${project_name}/list/${rule_name}`)).then(
      id_list => setIds(id_list)
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

    Object.assign(token, { seq, type, height, from: pos, to: pos + seq.length, filled });

    return [ pos + seq.length, height, new_filled ]
  }, [values]);

  const [ _, height ] = tokens.reduce( updateToken, [ 0, 0, true ] );

  const seq = tokens.map( t => t.seq ).join('');

  console.log({tokens})
  console.log({values})
  console.log({ids})

  const currentOptionsRegex = useMemo( () => {
    const option_sets = (tokens == null || tokens == []) ? null : tokens.map( (val, k) => {
      if (["resolved", "hidden"].includes(val.type)) {
        return val.seq
      } else if (values[k] != '') {
        return values[k]
      } else if (val.type == 'counter') {
        return '[0-9]+'
      } else {
        return '(' + Object.keys(val.values).join('|') + ')'
      }
    })
    if (option_sets != null) {
      return(''.concat(...option_sets))
    }
  }, [tokens, values])
  console.log({currentOptionsRegex})
  
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
      <IdsShower ids={matchIds(ids, currentOptionsRegex)}/>
  </Grid>
}

export default ComposeIdentifier;
