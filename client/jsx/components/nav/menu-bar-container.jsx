import * as ReactRedux from 'react-redux';
import MenuBar from './menu-bar';

const mapStateToProps = (state, ownProps)=>{

  // state == redux store
  return {

    appState: state['appState']
  };
}

const mapDispatchToProps = (dispatch, ownProps)=>{

  return {

    logOut: ()=>{

      var action = { type: 'LOG_OUT' };
      dispatch(action);
    }
  };
}

const MenuBarContainer = ReactRedux.connect(

  mapStateToProps,
  mapDispatchToProps,
)(MenuBar);

export default MenuBarContainer;