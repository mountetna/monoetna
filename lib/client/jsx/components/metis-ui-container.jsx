import * as ReactRedux from 'react-redux';
import MetisUI from './metis-ui';

const mapStateToProps = (state, ownProps)=>{
  // state == redux store
  return {
    user: state.user,
    files: state.files
  };
}

const mapDispatchToProps = (dispatch, ownProps)=>{
  return {};
}

const MetisUIContainer = ReactRedux.connect(
  mapStateToProps,
  mapDispatchToProps,
)(MetisUI);

export default MetisUIContainer;
