import * as ReactRedux from 'react-redux';
import LoginPanel from './login-panel';

const mapStateToProps = (state, ownProps)=>{

  // state == redux store
  return {

    appState: state['appState']
  };
}

const mapDispatchToProps = (dispatch, ownProps)=>{

  return {

    logIn: (email, pass)=>{

      var action = { type: 'LOG_IN', data: { email: email, pass: pass} };
      dispatch(action);
    }
  };
}

const LoginPanelContainer = ReactRedux.connect(

  mapStateToProps,
  mapDispatchToProps,
)(LoginPanel);

export default LoginPanelContainer;