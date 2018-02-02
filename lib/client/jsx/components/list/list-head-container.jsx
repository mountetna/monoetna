import * as ReactRedux from 'react-redux';
import ListHead from './list-head';

const mapStateToProps = (state, ownProps)=>{
  // state == redux store
  return {
    appState: state.appState
  };
}

const mapDispatchToProps = (dispatch, ownProps)=>{
  return {
    fileSelected: (fileObject)=>{
      let action = { type: 'FILE_SELECTED', fileObject };
      dispatch(action);
    }
  };
}

const ListHeadContainer = ReactRedux.connect(
  mapStateToProps,
  mapDispatchToProps,
)(ListHead);

export default ListHeadContainer;
