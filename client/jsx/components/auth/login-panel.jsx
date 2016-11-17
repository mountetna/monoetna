import * as React from 'react'

export default class LoginPanel extends React.Component{

  constructor(){

    super();
  }

  parseError(){

    if(this['props']['appState']['loginError']){

      return { display: 'block' };
    }
    else{

      return { display: 'none' };
    }
  }

  logIn(){

    var userEmail = document.getElementById('email-input').value;
    var userPass = document.getElementById('pass-input').value;
    this['props'].logIn(userEmail, userPass);
  }

  render(){

    return (

      <div id='login-group'>

        <input id='email-input' className='log-input' type='text' placeholder='Enter your email' />
        <br />
        <input id='pass-input' className='log-input' type='password' placeholder='Enter your password' />
        <br />
        <div className='log-error-message' style={ this.parseError() }>

          Invalid sign in.
        </div>
        <button className='login-button' onClick={ this['logIn'].bind(this) }>

          SIGN IN
        </button>
      </div>
    )
  }
}