import React, {useState, useEffect, useCallback, useRef} from 'react';
import Icon from 'etna-js/components/icon';
import {useReduxState} from 'etna-js/hooks/useReduxState';
import {selectUser} from 'etna-js/selectors/user-selector';
import { json_get, json_post } from 'etna-js/utils/fetch';
import { copyText } from 'etna-js/utils/copy';

import TokenBuilder from './token-builder.jsx';

import TextField from '@material-ui/core/TextField';
import Button from '@material-ui/core/Button';
import Select from '@material-ui/core/Select';
import MenuItem from '@material-ui/core/MenuItem';
import {makeStyles} from '@material-ui/core/styles';
import { isSuperuser } from 'etna-js/utils/janus';

const useStyles = makeStyles((theme) => ({
  text: {
    width: '300px'
  },
  keyText: {
    width: '650px'
  }
}));

const KeysSettings = ({user, setUser}) => {
  let [ error, setError ] = useState(null);
  const [disabled, setDisabled] = useState(false);

  const classes = useStyles();

  let { public_key } = user;

  let [ pem, setPemText ] = useState(null);

  let uploadKey = useCallback(
    () => {
      setDisabled(true);
      return json_post('/api/user/update_key', { pem })
      .then(
          ({user}) => { setUser(user); setError('Saved!'); }
      ).catch(
        e => e.then( ({error}) => setError(error) )
      ).finally(() => setDisabled(false));
    }, [pem]
  );

  return <div id='keys-group'>
    <div className='title'>Your Keys</div>
    {
      public_key
        ?
          <div className='item'>
            <i className='fa fa-key'></i> { public_key }
          </div>
        :
          <div className='item'>
            No registered keys
          </div>
    }
    <div className='item'>
      <TextField 
        multiline
        className={classes.keyText}
        onChange={ e => setPemText(e.target.value) }
        placeholder='Paste 2048+ bit RSA key in PEM format'/>
      <Button onClick={uploadKey} disabled={disabled}>Upload Key</Button>
      { error && <span className='error'>{error}</span> }
    </div>
  </div>;
};

const TaskTokenSettings = ({user}) => {
  let [ project_name, setProjectName ] = useState(null);

  let project_names = Object.keys(user.permissions).filter(
    project_name => project_name != 'administration'
  );

  let generateToken = useCallback(
    () => json_post(
      '/api/tokens/generate',
      { project_name, token_type: 'task' }
    ).then(
      (token) => copyText(token)
    ), [project_name]
  );
  const classes = useStyles();

  return <div id='task-token-group'>
    <div className='title'>Task Tokens</div>
    <div className='item'>
      <Select
        className={classes.text}
        value={ project_name }
        onChange={ e => setProjectName(e.target.value) }>
        {
          project_names.map(p => <MenuItem key={p} value={p}>{p}</MenuItem>)
        }
      </Select>
      <Button
        onClick={ project_name ? generateToken : null }
        disabled={ !project_name }>
        Copy Task Token
      </Button>
    </div>
  </div>;
};

const JanusSettings = () => {
  let user = useReduxState( state => selectUser(state) );

  let [ janusUser, setUser ] = useState({});

  useEffect(
    () => {
      json_get('/api/user/info').then(({user}) => setUser(user));
    }, []
  );

  return <div id='janus-settings'>
    <KeysSettings user={ janusUser } setUser={setUser} />
    <TaskTokenSettings user={ user }/>
    { isSuperuser(user) && <TokenBuilder user={ user }/> }
  </div>;
};

export default JanusSettings;
