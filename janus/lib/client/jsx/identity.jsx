import * as React from 'react';
import Cookies from 'js-cookie';
import { copyText } from 'etna-js/utils/copy';

import Button from '@material-ui/core/Button';

const copyToken = (e) => {
  let token = Cookies.get(CONFIG.token_name);
  if (token) copyText(token);
}

const Identity = ({user}) => 
  <div id='identity-group'>
    <div className='title'>Your Identity</div>
    <div className='item'>{ user.name }</div>
    <div className='item'>{ user.email }</div>
    <div className='item'><Button onClick={ copyToken }>Copy Token</Button></div>
  </div>;

export default Identity;
