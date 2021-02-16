import * as React from 'react';
import Cookies from 'js-cookie';

const copyToken = (e) => {
  let token = Cookies.get(CONFIG.token_name);
  if (token) {
    let input = document.createElement('input');
    input.value = token;
    document.body.appendChild(input);
    input.select();
    document.execCommand('copy',false);
    input.remove();
  }
}

const Identity = ({user}) => 
  <div id='identity-group'>
    <div className='title'>Your Identity</div>
    <div className='item'>{ user.name }</div>
    <div className='item'>{ user.email }</div>
    <div className='item'><button onClick={ copyToken }>Copy Token</button></div>
  </div>;

export default Identity;
