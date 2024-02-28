import React, {useState, useEffect, useCallback, useReducer} from 'react';
import { copyText, pasteText } from 'etna-js/utils/copy';
import { json_get, json_post } from 'etna-js/utils/fetch';
import { capitalize, jsonFormat } from 'etna-js/utils/format';
import { magmaPath } from 'etna-js/api/magma_api';

import Typography from '@material-ui/core/Typography';
import {makeStyles} from '@material-ui/core/styles';
import Grid from '@material-ui/core/Grid';
import Button from '@material-ui/core/Button';
import Paper from '@material-ui/core/Paper';
import TextField from '@material-ui/core/TextField';
import FormControl from '@material-ui/core/FormControl';
import Select from '@material-ui/core/Select';
import MenuItem from '@material-ui/core/MenuItem';
import Tooltip from '@material-ui/core/Tooltip';
import IconButton from '@material-ui/core/IconButton';
import HistoryIcon from '@material-ui/icons/HistoryRounded';
import AddIcon from '@material-ui/icons/Add';
import CodeIcon from '@material-ui/icons/Code';
import AutorenewIcon from '@material-ui/icons/Autorenew';

import Card from '@material-ui/core/Card';
import CardActionArea from '@material-ui/core/CardActionArea';
import CardActions from '@material-ui/core/CardActions';
import CardContent from '@material-ui/core/CardContent';

import ProjectHeader from 'etna-js/components/project-header';
import RevisionHistory from 'etna-js/components/revision-history';
import AddDialog from 'etna-js/components/add-dialog';

import RuleScript from './rule-script';

const upcase = v => v.replace(/[^A-Za-z0-9_]/g, '').toUpperCase();

const DEFAULT_STATE = {
  'rules': { },
  'synonyms': [],
  'tokens': {
    'SEP': {
      'name': 'SEP',
      'label': 'separator',
      'values': {
        '-': '# Separator'
      }
    },
  }
};

const removeRule = (rules, name) => {
  const { [name]: _, ...other_rules } = rules;

  return other_rules;
};

const removeToken = (tokens, name) => {
  const { [name]: _, ...other_tokens } = tokens;

  return other_tokens;
};

const removeSynFromSet = (state, action) => {
  const synonyms = state.synonyms || [];
  const update_set = [
    ...synonyms[action.pos].slice(0, action.ind),
    ...synonyms[action.pos].slice(action.ind+1)
  ];

  return [
    ...synonyms.slice(0,action.pos),
    ...update_set.length > 0 ? [ update_set ] : [],
    ...synonyms.slice(action.pos+1)
  ];
};

const addSynToSet = (state, action) => {
  const synonyms = state.synonyms || [];
  const update_set = [
    ...synonyms[action.pos].slice(0, action.ind),
    ...synonyms[action.pos].slice(action.ind+1)
  ];

  return [
    ...synonyms.slice(0,action.pos),
    [ ...synonyms[action.pos], action.token_name ],
    ...synonyms.slice(action.pos+1)
  ];
};

const reducer = (state, action) => {
  switch(action.type) {
    case 'SET':
      return action.state;
    case 'REMOVE_SYN_FROM_SET':
      return { ...state, synonyms: removeSynFromSet(state, action) };
    case 'ADD_SYN_TO_SET':
      return { ...state, synonyms: addSynToSet(state, action) };
    case 'ADD_SYNONYM_SET':
      return { ...state, synonyms: [ ...(state.synonyms || []), action.set ] };
    case 'ADD_RULE':
      return { ...state, rules: { ...state.rules, [action.name] : action.rule } };
    case 'ADD_TOKEN':
      return { ...state, tokens: { ...state.tokens, [action.token.name] : action.token } };
    case 'ADD_TOKEN_VALUE':
      return {
        ...state,
        tokens: {
          ...state.tokens,
          [action.name] : {
            ...state.tokens[action.name],
            values: {
              ...state.tokens[action.name].values,
              [ action.value ] : action.description
            }
          },
        }
      };
    case 'REMOVE_RULE':
      return { ...state, rules: removeRule(state.rules, action.name) };
    case 'REMOVE_TOKEN':
      return { ...state, tokens: removeToken(state.tokens, action.name) };
    default:
      return state;
  };
};

const useStyles = makeStyles((theme) => ({
  header: {
    borderBottom: '1px solid #ccc'
  },
  buttons: {
    '& > *': {
      margin: '0px 5px'
    },
    marginLeft: '40px',
    justifyContent: 'space-between',
    width: 'auto'
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
  rule_editor: {
    margin: '10px',
    border: '1px solid #ccc'
  },
  value: {
  },
  noprint: {
    '@media print': {
      display: 'none'
    }
  },
  strikeout: {
    '&:hover': {
      textDecoration: 'line-through',
      cursor: 'pointer',
      background: '#eee'
    }
  },
  loadingIcon: {
    'animation': '$spin 4s linear infinite'
  },
  '@keyframes spin': {
      '100%': {
          '-webkit-transform': 'rotate(360deg)',
          'transform': 'rotate(360deg)',
      }
  },
}));

const Synonym = ({set, pos, dispatch}) => {
  const classes = useStyles();

  return <Grid container className={classes.token}>
      {
        set.map(
          (name, ind) => <Grid className={classes.strikeout} onClick={ () => dispatch({ type: 'REMOVE_SYN_FROM_SET', pos, ind }) } key={name} item xs={1}>
            <Typography>{name}</Typography>
          </Grid>
        )
      }
      <Grid item xs={1}>
        <AddDialog
          buttonClass={classes.noprint}
          title='Add a synonym'
          content={`Add a new synonym for ${set.join(', ')}`}
          buttonText='ADD SYNONYM'
          update={ token_name => dispatch({ type: 'ADD_SYN_TO_SET', pos, token_name }) }
          mask={ upcase }
          placeholders={ [ 'TOKEN_NAME' ] }
        />
      </Grid>
  </Grid>;
};

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
    <Grid item xs={10}>{value}</Grid>
  </Grid>;
};

const Token = ({token, dispatch}) => {
  const classes = useStyles();

  return <Grid container className={classes.token}>
    <Grid item container direction='column' xs={2}>
      <Grid className={classes.strikeout} onClick={
        () => dispatch({type: 'REMOVE_TOKEN', name: token.name })
      } item>{token.name}</Grid>
      <Grid item><Typography color='secondary'>{token.label}</Typography></Grid>
    </Grid>
    <Grid item xs={10}>
      {
        Object.keys(token.values).map(
          name => <TokenValue key={name} name={name} value={token.values[name]} dispatch={dispatch} token={token}/>
        )
      }
      <AddDialog
        buttonClass={classes.noprint}
        title='Add a Value'
        content={
          <>
          </>
        }
        buttonText='ADD VALUE'
        update={ (value, description) => dispatch({ type: 'ADD_TOKEN_VALUE', name: token.name, value, description }) }
        mask={ v => v }
        placeholders={[ 'Value', 'Description' ]}
      />
    </Grid>
  </Grid>;
};

const Rule = ({name, rule, dispatch}) => {
  const classes = useStyles();

  return <Grid onClick={ () => dispatch({ type: 'REMOVE_RULE', name }) }
    container className={`${classes.rule} ${classes.strikeout}`}>
    <Grid item xs={2}>{ name }</Grid>
    <Grid item xs={10}>{ rule }</Grid>
  </Grid>;
};

const RuleEditorPane = ({type, update, mask, mask2, desc, placeholders, children}) => {
  const classes = useStyles();

  return <Card className={ classes.rule_editor } elevation={0}>
    <CardContent>
      <Typography gutterBottom variant='h5' component='h2'>
        {capitalize(type)+'s'}
      </Typography>
      { children }
    </CardContent>
    <CardActions>
      <AddDialog
        buttonClass={classes.noprint}
        title={`Add a ${capitalize(type)}`}
        content={desc}
        buttonText={ `ADD ${type.toUpperCase()}` }
        update={ update }
        mask={ mask }
        mask2={ mask2 }
        placeholders={ placeholders }
      />
    </CardActions>
  </Card>;
};

const RuleEditor = ({project_name}) => {
  const classes = useStyles();

  const [ editedState, dispatch ] = useReducer(reducer, DEFAULT_STATE);
  const [ savedState, setSavedState ] = useState(DEFAULT_STATE);
  const [ editedScript, setEditedScript ] = useState(jsonFormat(DEFAULT_STATE));
  const [ savedScript, setSavedScript ] = useState(jsonFormat(DEFAULT_STATE));
  const [ showRevisions, setShowRevisions ] = useState(null);
  const [ comment, setComment ] = useState('');
  const [ error, setError ] = useState('');
  const [ showJson, setShowJson ] = useState(false);

  const [loading, setLoading] = useState(true)

  const changed = showJson ? editedScript != savedScript : editedState != savedState;

  const setEditedState = state => dispatch({ type: 'SET', state });

  const unifyState = config => {
    const script = jsonFormat(config);
    setEditedState(config);
    setEditedScript(script);
    setSavedState(config);
    setSavedScript(script);
  };

  useEffect( () => {
    setLoading(true);
    json_get(magmaPath(`gnomon/${project_name}/`)).then(
      ({config}) => unifyState(config)
    );
    setLoading(false);
  }, [] );

  const saveRules = useCallback(
    () => {
      const newState = showJson ? JSON.parse(editedScript) : editedState;
      json_post(magmaPath(`gnomon/${project_name}`), { config: newState, comment }).then(
        () => {
          unifyState(newState);
          setError('');
        }
      ).catch(
        e => e.then( ({errors}) =>  setError( errors.join('; ')))
      );
    }, [ editedState, editedScript, savedState, showJson, comment ]
  );

  const revertRules = useCallback(
    () => {
      setEditedState(savedState);
      setEditedScript(savedScript);
      setError('');
      setComment('');

    }, [ editedState, savedState ]
  );

  return <Grid>
    <ProjectHeader project_name={project_name} className={classes.header}>
      <Grid container className={classes.buttons}>
      {
        changed && <>
            <Grid item>
              <TextField
                size='small'
                style={{width: 300}}
                value={comment}
                onChange={(e) => setComment(e.target.value)}
                placeholder='Revision comment'
              />
            </Grid>
            <Grid item>
              <Button
                disabled={comment == ''}
                onClick={ saveRules }>Save</Button>
            </Grid>
            <Grid item>
              <Button color='secondary' onClick={ revertRules }>Revert</Button>
            </Grid>
        </>
      }
      <Tooltip title={showJson ? 'show form' : 'show json'}>
        <IconButton
          disabled={changed}
          onClick={() => {
            setShowJson(!showJson);
            setError('');
          }}
          size='small'
          aria-label='revision history'
          color={showJson ? 'primary' : 'default'}
        >
          <CodeIcon />
        </IconButton>
      </Tooltip>
      <Tooltip title='revision history'>
        <IconButton
          className={classes.noprint}
          onClick={() => setShowRevisions(true)}
          size='small'
          aria-label='revision history'
        >
          <HistoryIcon />
        </IconButton>
      </Tooltip>
      {showRevisions != null && (
        <RevisionHistory
          getRevisions={() =>
            json_get(magmaPath(`gnomon/${project_name}/revisions`))
          }
          dateField='created_at'
          open={showRevisions}
          revisionDoc={(revision) =>
            JSON.stringify(revision.config, null, 2)
          }
          update={({config}) => {
            dispatch({ type: 'SET', state: config});
            setShowRevisions(false);
          }}
          onClose={() => setShowRevisions(false)}
        />
      )}
      {
        error && <Typography color='error'>{error}</Typography>
      }
      </Grid>
    </ProjectHeader>

    {
      loading ? <>
        <AutorenewIcon className={classes.loadingIcon}/>
        Loading
      </> : 
      showJson ? (<RuleScript script={editedScript} update={setEditedScript}/>) :
      <>
        <RuleEditorPane type='token'
          desc='Enter a token name and a label.'
          update={ (name, label) => dispatch({ type: 'ADD_TOKEN', token: { name, label, values: {} } }) }
          mask={ upcase }
          placeholders={['TOKEN_NAME', 'Token label']}>
          <Grid container className={classes.header}>
            <Grid item xs={2}>Name</Grid>
            <Grid item xs={10} container style={{ width: 'auto' }}>
              <Grid item xs={2}>Value</Grid>
              <Grid item xs={2}>Description</Grid>
            </Grid>
          </Grid>
          {
            Object.keys(editedState.tokens).map(
              name => <Token key={name} token={editedState.tokens[name]} dispatch={dispatch}/>
            )
          }
        </RuleEditorPane>
        <RuleEditorPane type='synonym'
          desc='Add a pair of synonyms.'
          update={ (name1,name2) => dispatch({ type: 'ADD_SYNONYM_SET', set: [ name1, name2 ]}) }
          mask={ upcase }
          mask2={ upcase }
          placeholders={['TOKEN_NAME1', 'TOKEN_NAME2']}>
          {
            editedState.synonyms?.map(
              (synonym_set,i) => <Synonym key={i} pos={i} set={synonym_set} dispatch={dispatch}/>
            )
          }
        </RuleEditorPane>
        <RuleEditorPane
          type='rule'
          desc={
            <>
              A <b>rule</b> is a sequence of space-separated token or rule names.
              <br/>
              Enter a name and a new rule.
            </>
          }
          update={ (name, rule) => dispatch({ type: 'ADD_RULE', name, rule }) }
          mask={ v => v.replace(/[^A-Za-z0-9_]/g, '').toLowerCase() }
          placeholders={['rule_name', 'Rule']}>
            <Grid container className={classes.header}>
              <Grid item xs={2}>Label</Grid>
              <Grid item xs={10}>Rule</Grid>
            </Grid>
            {
              Object.keys(editedState.rules).map(
                name => <Rule key={name} name={name} rule={editedState.rules[name]} dispatch={dispatch}/>
              )
            }
        </RuleEditorPane>
      </>
    }
  </Grid>;
};

export default RuleEditor;
