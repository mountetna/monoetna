import React, {useState, useEffect, useCallback, useReducer} from 'react';
import { copyText, pasteText } from 'etna-js/utils/copy';

import Typography from '@material-ui/core/Typography';
import {makeStyles} from '@material-ui/core/styles';
import Grid from '@material-ui/core/Grid';
import Button from '@material-ui/core/Button';
import TextField from '@material-ui/core/TextField';
import FormControl from '@material-ui/core/FormControl';
import Select from '@material-ui/core/Select';
import MenuItem from '@material-ui/core/MenuItem';
import IconButton from '@material-ui/core/IconButton';
import AddIcon from '@material-ui/icons/Add';

import Card from '@material-ui/core/Card';
import CardActionArea from '@material-ui/core/CardActionArea';
import CardActions from '@material-ui/core/CardActions';
import CardContent from '@material-ui/core/CardContent';

import Dialog from '@material-ui/core/Dialog';
import DialogActions from '@material-ui/core/DialogActions';
import DialogContent from '@material-ui/core/DialogContent';
import DialogContentText from '@material-ui/core/DialogContentText';
import DialogTitle from '@material-ui/core/DialogTitle';


import ProjectHeader from 'etna-js/components/project-header';

const DEFAULT_STATE = {
  "rules": {
    "project": "PROJECT",
    "subject": "PROJ SEP SP .n",
    "timepoint": ".subject SEP TP .n",
    "biospecimen": ".timepoint BSP .n",
    "immunoassay": ".biospecimen SEP IMM .n"
  },
  "tokens": {
    "PROJ": {
      "name": "PROJ",
      "label": "project",
      "values": {
        "MVIR1": "COMET"
      }
    },
    "PROJECT": {
      "name": "PROJECT",
      "label": "project",
      "values": {
        "mvir1": "COMET"
      }
    },
    "SP": {
      "name": "SP",
      "label": "species",
      "values": {
        "HS": "Homo sapiens"
      }
    },
    "BSP": {
      "name": "BSP",
      "label": "biospecimen",
      "values": {
        "PL": "Plasma",
        "SR": "Serum",
        "ETA": "Endotracheal aspirate",
        "BLD": "Whole blood",
        "NAS": "Nasal swab",
        "PBMC": "PBMCs"
      }
    },
    "IMM": {
      "name": "IMM",
      "label": "immunoassay",
      "values": {
        "CTK": "Cytokine",
        "VAG": "Viral Antigen",
        "LNK": "Olink"
      }
    },
    "TP": {
      "name": "TP",
      "label": "timepoint",
      "values": {
        "D": "Day",
        "DN": "Negative Day",
        "M": "Month"
      }
    },
    "SEP": {
      "name": "SEP",
      "label": "separator",
      "values": {
        "-": "# Separator"
      }
    }
  }
};

const reducer = (state, action) => {
  switch(action.type) {
    case 'PASTE':
      return {
        rules: { ...action.paste.rules },
        tokens: { ...action.paste.tokens }
      };
    case 'ADD_RULE':
      return { ...state, rules: { ...state.rules, [action.name] : action.rule } }
    case 'ADD_TOKEN':
      return { ...state, tokens: { ...state.tokens, [action.token.name] : action.token } }
    case 'REMOVE_RULE':
      const { [action.name]: _, ...other_rules } = state.rules;
      return { ...state, rules: other_rules };
    default:
      return state;
  }
};

const AddToken = ({dispatch}) => {
  const [ open, setOpen ] = useState(false);
  const [ name, setName ] = useState('');
  const [ label, setLabel ] = useState('');

  const handleClose = () => {
    setOpen(false);
    setName('');
    setLabel('');
  }

  const handleAdd = useCallback( () => {
    dispatch({ type: 'ADD_TOKEN', token: { name, label, values: {} } });
    handleClose();
  }, [name, label]);

  return <React.Fragment>
    <Button 
    variant='text'
    startIcon={<AddIcon/>}
    title='Add Token'
    color='secondary'
    onClick={() => setOpen(true)}>ADD TOKEN</Button>
    <Dialog open={open} onClose={handleClose} aria-labelledby='form-dialog-title'>
      <DialogTitle id='form-dialog-title'>Add a Token</DialogTitle>
      <DialogContent>
        <DialogContentText>Enter a token name and a label.</DialogContentText>
        <TextField
          autoFocus
          margin='dense'
          placeholder='TOKEN_NAME'
          fullWidth
          value={ name }
          onChange={ e => setName(e.target.value.replace(/[^A-Za-z0-9_]/g, '').toUpperCase()) }
        />
        <TextField
          margin='dense'
          placeholder='Token label'
          fullWidth
          value={ label }
          onChange={ e => setLabel(e.target.value) }
        />
      </DialogContent>
      <DialogActions>
        <Button onClick={handleClose} color='primary'>
          Cancel
        </Button>
        <Button onClick={ handleAdd } color='primary'>
          Add
        </Button>
      </DialogActions>
    </Dialog>
  </React.Fragment>
}

const useStyles = makeStyles((theme) => ({
  header: {
    borderBottom: '1px solid #eee'
  },
  token: {
    borderBottom: '1px solid #eee',
    padding: '10px 0px'
  },
  rule: {
    borderBottom: '1px solid #eee',
    padding: '10px 0px'
  },
  row: {
    flex: '60px'
  },
  rows: {
    flex: '1 1 auto'
  },
  editor: {
    marginLeft: 20,
    height: 'calc(100vh - 61px - 48px)',
    position: 'relative',
    overflowY: 'clip'
  },
  value: {
  },
  strikeout: {
    '&:hover': {
      textDecoration: 'line-through',
      cursor: 'pointer',
      background: '#eee'
    }
  }
}));


const TokenValue = ({name,value,dispatch, token}) => {
  const classes = useStyles();
  return <Grid
    onClick={
      () => {
        const { [name]: _, ...other_values } = token.values;
        dispatch({
          type: 'ADD_TOKEN',
          token: {
            ...token,
            values: other_values
          }
        });
      }
    }
    className={`${classes.value} ${classes.strikeout}`} container>
    <Grid item xs={2}>{name}</Grid>
    <Grid item xs={2}>{value}</Grid>
  </Grid>
}

const AddTokenValue = ({token, dispatch}) => {
  const [ open, setOpen ] = useState(false);
  const [ value, setValue ] = useState('');
  const [ description, setDescription ] = useState('');

  const handleClose = () => {
    setOpen(false);
    setValue('');
    setDescription('');
  }

  const handleAdd = useCallback( () => {
    dispatch({ type: 'ADD_TOKEN', token: {
      ...token, values: {
        ...token.values,
        [ value ] : description
      }
    }});
    handleClose();
  }, [ value, description]);

  return <React.Fragment>
    <Button 
    variant='text'
    startIcon={<AddIcon/>}
    title='Add Token'
    color='secondary'
    onClick={() => setOpen(true)}>ADD VALUE</Button>
    <Dialog open={open} onClose={handleClose} aria-labelledby='form-dialog-title'>
      <DialogTitle id='form-dialog-title'>Add a Value</DialogTitle>
      <DialogContent>
        <DialogContentText>Enter a possible value and description for token {token.name}.</DialogContentText>
        <TextField
          autoFocus
          margin='dense'
          placeholder='Value'
          fullWidth
          value={ value }
          onChange={ e => setValue(e.target.value) }
        />
        <TextField
          margin='dense'
          placeholder='Description'
          fullWidth
          value={ description }
          onChange={ e => setDescription(e.target.value) }
        />
      </DialogContent>
      <DialogActions>
        <Button onClick={handleClose} color='primary'>
          Cancel
        </Button>
        <Button onClick={ handleAdd } color='primary'>
          Add
        </Button>
      </DialogActions>
    </Dialog>
  </React.Fragment>
}

const Token = ({token, dispatch}) => {
  const classes = useStyles();

  return <Grid container className={classes.token}>
    <Grid item container direction='column' xs={2}>
      <Grid item>{token.name}</Grid>
      <Grid item><Typography color='secondary'>{token.label}</Typography></Grid>
    </Grid>
    <Grid item xs={10}>
      {
        Object.keys(token.values).map(
          name => <TokenValue key={name} name={name} value={token.values[name]} dispatch={dispatch} token={token}/>
        )
      }
      <AddTokenValue token={token} dispatch={dispatch}/>
    </Grid>
  </Grid>
}

const AddRule = ({dispatch}) => {
  const [ open, setOpen ] = useState(false);
  const [ name, setName ] = useState('');
  const [ rule, setRule ] = useState('');

  const handleClose = () => {
    setOpen(false);
    setName('');
    setRule('');
  }

  const handleAdd = useCallback( () => {
    dispatch({ type: 'ADD_RULE', name, rule });
    handleClose();
  }, [name, rule]);

  return <React.Fragment>
    <Button 
    variant='text'
    startIcon={<AddIcon/>}
    title='Add Rule'
    color='secondary'
    onClick={() => setOpen(true)}>ADD RULE</Button>
    <Dialog open={open} onClose={handleClose} aria-labelledby='form-dialog-title'>
      <DialogTitle id='form-dialog-title'>Add a Rule</DialogTitle>
      <DialogContent>
        <DialogContentText>Enter a name and a new rule. A <b>rule</b> is a sequence of space-separated token or rule names.</DialogContentText>
        <TextField
          autoFocus
          margin='dense'
          placeholder='rule_name'
          fullWidth
          value={ name }
          onChange={ e => setName(e.target.value.replace(/[^A-Za-z0-9_]/g, '').toLowerCase()) }
        />
        <TextField
          margin='dense'
          placeholder='Rule'
          fullWidth
          value={ rule }
          onChange={ e => setRule(e.target.value) }
        />
      </DialogContent>
      <DialogActions>
        <Button onClick={handleClose} color='primary'>
          Cancel
        </Button>
        <Button onClick={ handleAdd } color='primary'>
          Add
        </Button>
      </DialogActions>
    </Dialog>
  </React.Fragment>
}

const Rule = ({name, rule, dispatch}) => {
  const classes = useStyles();

  return <Grid onClick={ () => dispatch({ type: 'REMOVE_RULE', name }) }
    container className={`${classes.rule} ${classes.strikeout}`}>
    <Grid item xs={2}>{ name }</Grid>
    <Grid item xs={10}>{ rule }</Grid>
  </Grid>
}

const RuleEditor = ({project_name}) => {
  const classes = useStyles();

  const [ state, dispatch ] = useReducer(reducer, DEFAULT_STATE);

  return <Grid>
    <ProjectHeader project_name={project_name} className={classes.header}/>
    <Card elevation={0}>
      <CardContent>
        <Typography gutterBottom variant="h5" component="h2">
          Tokens
        </Typography>
        <Grid container className={classes.header}>
          <Grid item xs={2}>Name</Grid>
          <Grid item xs={10} container style={{ width: 'auto' }}>
            <Grid item xs={2}>Value</Grid>
            <Grid item xs={2}>Description</Grid>
          </Grid>
        </Grid>
        {
          Object.keys(state.tokens).map(
            name => <Token key={name} token={state.tokens[name]} dispatch={dispatch}/>
          )
        }
      </CardContent>
      <CardActions>
        <AddToken dispatch={dispatch}/>
      </CardActions>
    </Card>
    <Card elevation={0}>
      <CardContent>
        <Typography gutterBottom variant="h5" component="h2">
          Rules
        </Typography>
        <Grid container className={classes.header}>
          <Grid item xs={2}>Label</Grid>
          <Grid item xs={10}>Rule</Grid>
        </Grid>
        {
          Object.keys(state.rules).map(
            name => <Rule key={name} name={name} rule={state.rules[name]} dispatch={dispatch}/>
          )
        }
      </CardContent>
      <CardActions>
        <AddRule dispatch={dispatch}/>
      </CardActions>
    </Card>
    <Button onClick={
      () => copyText( JSON.stringify(state,null,2) )
    }>Copy</Button>
    <TextField inputProps={{ id: 'paste' }}/>
    <Button onClick={
      () => {
        const paste = JSON.parse(document.querySelector('#paste').value);
        dispatch({ type: 'PASTE', paste })
      }
    }>Paste</Button>
  </Grid>
}

export default RuleEditor;
