import * as React from 'react';
import Cookies from 'js-cookie';
import { copyText } from 'etna-js/utils/copy';

const copyToken = (e) => {
  let token = Cookies.get(CONFIG.token_name);
  if (token) copyText(token);
}

const Identity = ({user}) => 
  <div id='identity-group'>
    <div className='title'>Your Identity</div>
    <div className='item'>{ user.name }</div>
    <div className='item'>{ user.email }</div>
    <div className='item'><button onClick={ copyToken }>Copy Token</button></div>
  </div>;

export default Identity;
