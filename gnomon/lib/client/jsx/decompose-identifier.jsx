import React from 'react';
import Grid from '@material-ui/core/Grid';
import Typography from '@material-ui/core/Typography';
import ProjectHeader from 'etna-js/components/project-header';
import {makeStyles} from '@material-ui/core/styles';

const useStyles = makeStyles((theme) => ({
  header: {
    borderBottom: '1px solid #eee'
  },
  tokens: {
    width: '100%',
    position: 'absolute'
  },
  letters: {
    width: '100%'
  },
  rules: {
    width: '100%',
    position: 'absolute'
  },
  decomposer: {
    marginLeft: 20,
    height: 'calc(100vh - 61px - 48px)',
    position: 'relative',
    overflowY: 'clip'
  },
  letter_box: {
    position: 'relative',
    width: '40px',
    height: '40px',
    border: '1px solid #ccc',
    background: '#eee',
    borderRight: 'none',
    '$counter + &:not($counter)': {
      borderLeft: '1px solid #888'
    }
  },
  separator_box: {
    position: 'relative',
    width: '40px',
    height: '40px',
    padding: '2px',
    border: '1px solid transparent',
    '$counter + &': {
      borderLeft: '1px solid #888',
    },
    ':not($counter) + &': {
      borderLeft: '1px solid #ccc',
    }
  },
  letter: {
    fontSize: '30px'
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

const isSeparator = letter => letter.match(/[\.\-\_]/)

const Rule = ({rule, order, total, identifier}) => {
  const voff = 50;
  return <Grid style={{ position: 'absolute', left: 0, top: 0 }}>
    {/* bracket */}
    <Grid style={{ position: 'absolute',
      top: voff + 2.5 + order * 25,
      left: rule.from * 40,
      width: (rule.to - rule.from + 1) * 40, border: '1px solid black',
      height: '10px', borderTop: 'none' }}/>
    {/* text */}
    <Grid style={{
      whiteSpace: 'nowrap',
      position: 'absolute',
      top: voff - 5 + 25 * order,
      left: 50 + identifier.length * 40 }}>
      { rule.label }
    </Grid>
    {/* pointer */}
    <Grid style={{
      position: 'absolute',
      opacity: 0.1,
      width: identifier.length * 40 - (rule.to) * 40 + 5,
      left: (1+rule.to) * 40 + 1,
      top: voff + 7.5 + 25 * order,
      borderTop: '1px solid black'
      }}/>
  </Grid>
}

const Token = ({token, order, total, identifier}) => {
  const voff = 40;
  return <Grid style={{ position: 'absolute', left: 0, top: 0 }}>
    {/* bracket */}
    <Grid style={{
      position: 'absolute',
      bottom: voff,
      left: token.from * 40,
      width: (token.to - token.from + 1) * 40, border: '1px solid black',
      height: '10px', borderBottom: 'none' }}/>
    {/* text */}
    <Grid style={{
      whiteSpace: 'nowrap',
      position: 'absolute',
      bottom: voff + 25 * (total - order),
      left: 50 + identifier.length * 40 }}>
      { token.label } : { token.description }
    </Grid>
    {/* pointer */}
    <Grid style={{
      position: 'absolute',
      width: identifier.length * 40 - (token.to + token.from + 1) * 20 + 45,
      left: (token.to + token.from + 1) * 20,
      bottom: voff + 10, borderLeft: '1px solid black',
      opacity: 0.1,
      borderTop: '1px solid black',
      height: 25 * (total - order)  }}/>
  </Grid>
}

const Letter = ({letter, position, counter}) => {
  const classes = useStyles();
  return <Grid className={
      `${isSeparator(letter) ? `sep ${classes.separator_box}` : classes.letter_box} ${counter ? `count ${classes.counter}` : ''}`
    }
    container alignItems='center' justify='center'>
    <span className={classes.letter}>{letter}</span>
    <span className={classes.position}>{position}</span>
  </Grid>
}

const DecomposeIdentifier = ({project_name, identifier}) => {
  const classes = useStyles();

  identifier = 'MVIR1-HS169-D0PL1-CTK1';

  const decomposition = {
    rules: [
      { from: 0, to: 4, label: 'project' },
      { from: 0, to: 10, label: 'subject' },
      { from: 0, to: 13, label: 'timepoint' },
      { from: 0, to: 16, label: 'biospecimen' },
      { from: 0, to: 21, label: 'immunoassay' },
    ],
    tokens: [
      { from: 0, to: 4, label: 'project', description: 'COMET' },
      { from: 6, to: 7, label: 'species', description: 'Homo sapiens' },
      { from: 12, to: 12, label: 'day', description: 'Day' },
      { from: 14, to: 15, label: 'biospecimen_type', description: 'Plasma' },
      { from: 18, to: 20, label: 'assay_type', description: 'Cytokine' }
    ],
    counters: '0000000011100100100001'
  }

  return <Grid>
    <ProjectHeader project_name={ project_name } className={classes.header}/>
    <Grid container alignItems='center' justify='center' className={classes.decomposer} style={{ width: 40 * identifier.length }}>
      <Grid container alignItems='center' justify='center' className={classes.tokens}>
        {
          decomposition.tokens.map(
            (token,i) => <Token key={i} token={token} order={i} total={ decomposition.tokens.length } identifier={identifier}/>
          )
        }
      </Grid>
      {
        identifier.split('').map(
          (l,i) => <Letter key={i} position={i} letter={l} counter={ decomposition.counters[i] == '1' } />
        )
      }
      <Grid container alignItems='center' justify='center' className={classes.rules}>
        {
          decomposition.rules.map(
            (rule,i) => <Rule key={i} rule={rule} order={i} total={ decomposition.rules.length } identifier={identifier}/>
          )
        }
      </Grid>
    </Grid>
  </Grid>
}

export default DecomposeIdentifier;
