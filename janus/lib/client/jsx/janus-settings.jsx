import React, {useState, useEffect, useCallback, useRef} from 'react';
import Icon from 'etna-js/components/icon';
import SelectInput from 'etna-js/components/inputs/select_input';
import {useReduxState} from 'etna-js/hooks/useReduxState';
import {selectUser} from 'etna-js/selectors/user-selector';
import { json_get, json_post } from 'etna-js/utils/fetch';
import { copyText } from 'etna-js/utils/copy';

const KeysSettings = ({user}) => {
  let [ error, setError ] = useState(null);

  let { public_key } = user;

  let pemText = useRef(null);

  let uploadKey = useCallback(
    () => json_post('/update_key', { pem: pemText.current.value }).then(
      ({user}) => { setUser(user); setError(null); }
    ).catch(
      e => e.then( ({error}) => setError(error) )
    ), [pemText]
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
      <textarea ref={pemText}
        placeholder='Paste 2048+ bit RSA key in PEM format'></textarea>
      <button onClick={uploadKey}>Upload Key</button>
      { error && <span className='error'>{error}</span> }
    </div>
  </div>
}

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

  return <div id='task-token-group'>
    <div className='title'>Task Tokens</div>
    <div className='item'>
      <SelectInput
        values={ project_names }
        value={ project_name }
        onChange={ setProjectName }
        showNone='enabled'
      />
      <button
        onClick={ project_name ? generateToken : null }
        disabled={ !project_name }>
        Copy Task Token
      </button>
    </div>
  </div>;
}

const JanusSettings = () => {
  let user = useReduxState( state => selectUser(state) );

  let [ janusUser, setUser ] = useState({});

  useEffect(
    () => {
      json_get('/user').then(({user}) => setUser(user))
    }, []
  );

  return <div id='janus-settings'>
    <KeysSettings user={ janusUser }/>
    <TaskTokenSettings user={ user }/>
  </div>;
}

export default JanusSettings;
