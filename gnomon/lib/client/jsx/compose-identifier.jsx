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
import Letter from './letter';
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
  }
}));

const TokenEditor = ({token, seq, height, value, pos, update}) => {
  const classes = tokenEditorStyles()
  if (token.type == 'hidden') return null;
  const contents =(token.type == 'resolved')
  ?  <TextField label={token.label} InputLabelProps={{ shrink: true }}  InputProps={{ readOnly: true }} value={ firstValue(token) }/>
  : (token.type == 'counter')
  ? <><TextField onChange={ e => update(pos, e.target.value) } value={value} InputLabelProps={{ shrink: true }} size='small' inputProps={{ pattern: "[0-9]*" }} label={token.label}/> <IconButton title={ `Fill the next available value for ${token.label}`}><PlusOneIcon/></IconButton></>
  : (token.type == 'unresolved')
  ? <FormControl style={{ width: 230 }}>
      <InputLabel shrink>{token.label}</InputLabel>
      <Select onChange={ e => update(pos, e.target.value) } value={value} displayEmpty size='small'>
        <MenuItem value=''><em>None</em></MenuItem>
        {
          Object.keys(token.values).map(
            name => <MenuItem key={name} value={name}>{token.values[name]}</MenuItem>
          )
        }
      </Select>
    </FormControl>
  : null

  const voff = 25;
  const lh = 70;

  const h = voff + 20 + (15 - lh / 2) * height + (lh - 15) * token.height;
  const w = 75;
  const r = Math.sqrt(h*h + w*w);
  const a = Math.atan(h/w);


  return <Grid className={classes.token_editor} style={{position: 'absolute', left: 0, top: 0}}>
    {/* bracket */}
    <Grid style={{ position: 'absolute',
      bottom: voff + 2.5,
      left: token.from * 40+1,
      width: (token.to - token.from + 1) * 40 - 4, border: '1px solid black',
      height: '10px', borderBottom: 'none' }}/>
    {/* text */}
    <Grid container alignItems='center' style={{
      width: 290,
      position: 'absolute',
      left: w + 25 + seq.length * 40,
      bottom: (height / 2 - token.height) * lh - 2,
      border: '1px solid #888',
      padding: '5px',
      whiteSpace: 'nowrap',
      background: token.type == 'resolved' ? 'none' : '#eee'
    }}>
    {contents}
    </Grid>
    {/* line1 */}
    <Grid
      className={ classes.pointer }
      style={{
      height: 0,
      left: (token.to + token.from + 1) * 20,
      bottom: voff + 12.5,
      transformOrigin: 'bottom left',
      transform: 'rotate(-90deg)',
      width: 15 * (height - token.height + 1)  }}/>
    {/* line2 */}
    <Grid
      className={ classes.pointer }
      style={{
      width: seq.length * 40 - (token.to + token.from + 1) * 20,
      left: (token.to + token.from + 1) * 20,
      bottom: voff + 12 + 15 * (height - token.height + 1),
      height: 0 }}/>
    {/* line3 */}
    <Grid 
      className={ classes.pointer }
      style={{
        transformOrigin: 'bottom left',
        left: seq.length * 40,
        width: r,
        transform: `rotate(${a}rad)`,
        bottom: voff + 12 + 15 * (height - token.height + 1),
        height: 0  }}/>
    {/* line4 */}
    <Grid
      className={ classes.pointer }
      style={{
        left: seq.length * 40 + w,
        width: 20,
        bottom: 7 + (height / 2 - token.height) * lh,
        height: 0  }}/>
  </Grid>
}

const Token = ({token}) => {
  const classes = useStyles();

  return <Grid style={{position: 'relative'}}>
    <Letters className={classes[token.type]} seq={ token.seq }/>
  </Grid>
}

const ComposeIdentifier = ({project_name}) => {
  const classes = useStyles();

  // 'MVIR1-HS169-D0PL1-CTK1';

  // a string of tokens we must satisfy
  const model = {
    tokens: [
      {
        name: 'PROJ',
        label: 'project',
        values: { 'MVIR1': 'COMET' }
      },
      {
        name: 'SEP',
        values: { '-' : '# Separator' }
      },
      {
        name: 'SP',
        label: 'species',
        values: { 'HS' : 'Homo Sapiens' }
      },
      {
        name: 'n',
        label: 'Subject counter'
      },
      {
        name: 'SEP',
        values: { '-' : '# Separator' }
      },
      {
        name: 'TP',
        label: 'timepoint',
        values: { 'D' : 'Day', 'DN' : 'Negative Day', 'M' : 'Month' }
      },
      {
        name: 'n',
        label: 'Timepoint counter'
      },
      {
        name: 'BSP',
        label: 'biospecimen',
        values: {
          'PL' : 'Plasma',
          'BLD' : 'Blood',
          'PBMC' : 'PBMCs',
          'SR' : 'Serum',
          'ETA' : 'Endotracheal Aspirate'
        }
      },
      {
        name: 'n',
        label: 'Biospecimen counter'
      },
      {
        name: 'SEP',
        values: { '-' : '# Separator' }
      },
      {
        name: 'IMM',
        label: 'immunoassay',
        values: {
          'CTK': 'Cytokine',
          'VAG' : 'Viral Antigen',
          'PDV' : 'Viral PhIP-Seq',
          'LNK' : 'Olink',
          'RSL' : 'Recan Luminex'
        }
      },
      {
        name: 'n',
        label: 'Immunoassay counter'
      },
    ]
  };

  const [ values, setValues ] = useState(model.tokens.map(t=>''));

  const setValue = useCallback( (i, val) => {
    let newValues = [ ... values ];
    newValues[i] = val;
    console.log({newValues});
    setValues( newValues );
  }, [ values ]);

  const updateToken = useCallback(([ pos, height ], token, i) => {
    const [ seq, type ] = (token.values && Object.keys(token.values).length == 1)
      ? (firstValue(token)[0] == '#'
         ? [ firstKey(token), 'hidden' ]
         : [ firstKey(token), 'resolved' ] )
      : token.name == 'n'
      ? [ values[i] || token.name, 'counter' ]
      : [ values[i] || token.name, 'unresolved' ];

    if (type != 'hidden') height = height + 1;

    Object.assign(token, { seq, type, height, from: pos, to: pos + seq.length - 1 });

    return [ pos + seq.length, height ]
  }, [values]);

  const [ _, height ] = model.tokens.reduce( updateToken, [ 0, 0 ] );

  const seq = model.tokens.map( t => t.seq ).join('');

  return <Grid>
    <ProjectHeader project_name={ project_name } className={classes.header}/>
    <Grid container alignItems='center' className={classes.composer} style={{ width: 40 * (seq.length+1) }}>
      <Grid container className={classes.tokens}>
      {
        model.tokens.map(
          (token, i) => <TokenEditor key={i} token={token} seq={seq} height={height} update={setValue} value={ values[i] } pos={i}/>
        )
      }
      </Grid>
      {
        model.tokens.map(
          (token,i) => <Token key={i} token={token}/>
        )
      }
    </Grid>
  </Grid>
}

export default ComposeIdentifier;
