import * as ReactRedux from 'react-redux';
import BasicView from './basic-view';

const mapStateToProps = (state, ownProps)=>{

  // state == redux store
  return {

    'userInfo': state['userInfo']
  };
}

const mapDispatchToProps = (dispatch, ownProps)=>{

  return {};
}

const BasicViewContainer = ReactRedux.connect(

  mapStateToProps,
  mapDispatchToProps,
)(BasicView);

export default BasicViewContainer;