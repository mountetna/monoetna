import React, {useState, useEffect } from 'react';
import Grid from '@material-ui/core/Grid';
import Link from '@material-ui/core/Link';
import Typography from '@material-ui/core/Typography';
import ProjectHeader from 'etna-js/components/project-header';
import {makeStyles} from '@material-ui/core/styles';
import Letter from './letter';
import Bracket from './bracket';
import Corner from './corner';
import { json_get, json_post } from 'etna-js/utils/fetch';
import { dateFormat } from 'etna-js/utils/format';
import { magmaPath } from 'etna-js/api/magma_api';
import {getDocuments} from 'etna-js/api/magma_api';
import TokenEditor from './token-editor';
import CheckBoxOutlinedIcon from '@material-ui/icons/CheckBoxOutlined';
import LaunchIcon from '@material-ui/icons/Launch';
import Table from '@material-ui/core/Table';
import TableBody from '@material-ui/core/TableBody';
import TableCell from '@material-ui/core/TableCell';
import TableContainer from '@material-ui/core/TableContainer';
import TableHead from '@material-ui/core/TableHead';
import TableRow from '@material-ui/core/TableRow';
import Paper from '@material-ui/core/Paper';
import Tooltip from '@material-ui/core/Tooltip';
import { LinkedIdTable } from './idTreeTable';

const useStyles = makeStyles((theme) => ({
  header: {
    borderBottom: '1px solid #eee'
  },
  tokens: {
    width: '100%',
    position: 'absolute'
  },
  letters: {
    width: '100%',
    position: 'absolute',
    top: -100
  },
  rules: {
    width: '100%',
    position: 'absolute',
    top: 0
  },
  decomposer: {
    marginLeft: 20,
    height: 'calc(100vh - 61px - 48px)',
    position: 'relative',
    overflowY: 'clip'
  },
  letter: {
    border: '1px solid #ccc',
    background: '#eee',
    borderRight: 'none',
    '$counter + &:not($counter)': {
      borderLeft: '1px solid #888'
    }
  },
  separator: {
    padding: '2px',
    border: '1px solid transparent',
    '$counter + &': {
      borderLeft: '1px solid #888',
    },
    ':not($counter) + &': {
      borderLeft: '1px solid #ccc',
    }
  },
  counter: {
    background: 'rgba(128,128,128,0.5)',
    border: '1px solid #888',
    '&:not(:last-of-type)': {
      borderRight: 'none'
    }
  },
  position: {
    position: 'absolute',
    top: '40px',
    width: '40px',
    textAlign: 'center',
    fontSize: '10px',
    color: '#444'
  },
  token: {
    position: 'absolute',
    bottom: '50px'
  }
}));

const isSeparator = token => token == 'SEP';
const isCounter = token => token.match(/_counter/);

const tokenLabel = (token_name, grammar) => (
  token_name in grammar.tokens ? grammar.tokens[token_name].label : token_name
);
const tokenDescription = (token_name, value, grammar) => (
  token_name in grammar.tokens ? grammar.tokens[token_name].values[value] : 'Numeric counter'
);

const DecomposeIdentifier = ({project_name, identifier}) => {
  const classes = useStyles();

  const [ decomposition, setDecomposition ] = useState(null);
  const [ grammar, setGrammar ] = useState(null);
  const [ models, setModels ] = useState([]);

  useEffect( () => {
    json_get(magmaPath(`gnomon/${project_name}/decompose/${identifier}`)).then(
      decomposition => setDecomposition(decomposition)
    );
    json_get(magmaPath(`gnomon/${project_name}/`)).then(
      ({config}) => {
        setGrammar(config);
      }
    );
    getDocuments(
      {
        project_name,
        model_name: 'all',
        record_names: [],
        attribute_names: 'all'
      },
      fetch
    )
      .then(({models}) => setModels(Object.keys(models)))
      .catch((e) => console.log({e}));
  }, [] );

  let from = 0;
  let to = 0;
  let height = -1;
  const total = decomposition ? decomposition.tokens.filter(([token,seq]) => !isSeparator(token)).length : 0;

  return <Grid>
    <ProjectHeader project_name={ project_name } className={classes.header}/>
    <Grid container alignItems='center' justify='center' className={classes.decomposer} style={{ width: 40 * identifier.length }}>
      <Grid container alignItems='center' justify='center' className={classes.tokens}>
        {
          grammar && decomposition && decomposition.tokens.map(
            ([token,seq],i) => {
              from = to;
              to = from + seq.length;

              if (isSeparator(token)) return <div/>;

              height = height + 1;

              return <TokenEditor
                key={i}
                token={{ from, to, height, type: 'resolved', label: tokenLabel(token, grammar)}}
                description={ tokenDescription(token, seq, grammar) }
                seq={identifier}
                height={total}
                pos={i}
                voff={120}/>
            }
          )
        }
        <Grid item container className={classes.letters}>
        {
          decomposition && decomposition.tokens.map(
            ([token, seq], i) => seq.split('').map( 
              (l,j) => <Letter key={`${i}.${j}`}
                className={
                  `${
                    isSeparator(token) ? classes.separator : classes.letter
                  } ${
                    isCounter(token) ? classes.counter : ''
                  }`
                }
                letter={l} />
            )
          ).flat()
        }
        </Grid>
          { decomposition && <LinkedIdTable decomposition={decomposition} project_name={project_name}/> }
      </Grid>
    </Grid>
  </Grid>
}

export default DecomposeIdentifier;
