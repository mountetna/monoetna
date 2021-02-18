import React, {useState, useEffect, useCallback, useRef} from 'react';
import Icon from 'etna-js/components/icon';
import {useReduxState} from 'etna-js/hooks/useReduxState';
import {selectUser} from 'etna-js/selectors/user-selector';
import { json_get, json_post } from 'etna-js/utils/fetch';

const JanusSettings = () => {
  let user = useReduxState( state => selectUser(state) );

  let [ error, setError ] = useState(null);

  let [ janusUser, setUser ] = useState({});

  useEffect(
    () => {
      json_get('/user').then(({user}) => setUser(user))
    }, []
  );

  let { public_key } = janusUser;

  let pemText = useRef(null);

  let uploadKey = useCallback(
    () => json_post('/update_key', { pem: pemText.current.value }).then(
      ({user}) => { setUser(user); setError(null); }
    ).catch(
      e => e.then( ({error}) => setError(error) )
    ), [pemText]
  );

  return <div id='janus-settings'>
    <div id='keys-group'>
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
  </div>;
}

export default JanusSettings;
