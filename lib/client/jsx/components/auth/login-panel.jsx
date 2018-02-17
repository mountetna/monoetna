import * as React from 'react'

export default class LoginPanel extends React.Component{
  constructor(){
    super();
  }

  parseError(){
    if(this.props.user.loginError){
      return { display: 'block' };
    }
    else{
      return { display: 'none' };
    }
  }

  logIn(){
    var userEmail = document.getElementById('email-input').value;

    if(!VALIDATE_EMAIL(userEmail)){
      this.props.user.loginError = true;
      this.props.user.loginErrorMsg = 'Bad email.'
      this.forceUpdate()
      return;
    }

    var userPass = document.getElementById('pass-input').value;
    this.props.logIn(userEmail, userPass);
  }

  runOnEnter(event){
    event = event || window.event;
    if(event.keyCode == 13 || event.which == 13){
      var userEmail = document.getElementById('email-input').value;
      var userPass = document.getElementById('pass-input').value;

      if(userPass != '' && VALIDATE_EMAIL(userEmail)){
        this.props.logIn(userEmail, userPass);
      }
    }
    else{
      return;
    }
  }

  render(){
    return (
      <div id='login-group'>
        <input id='email-input' className='log-input' type='text' placeholder='Enter your email' onKeyPress={ this.runOnEnter.bind(this) }/>
        <br />
        <input id='pass-input' className='log-input' type='password' placeholder='Enter your password' onKeyPress={ this.runOnEnter.bind(this) }/>
        <br />
        <div className='log-error-message' style={ this.parseError() }>
          { this.props.user.loginErrorMsg }
        </div>
        <button className='login-button' onClick={ this.logIn.bind(this) }>
          SIGN IN
        </button>
      </div>
    )
  }
}
