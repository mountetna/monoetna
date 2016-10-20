import * as ReactRedux from 'react-redux';
import ListHead from './list-head';

const mapStateToProps = (state, ownProps)=>{

  // state == redux store
  return {};
}

const mapDispatchToProps = (dispatch, ownProps)=>{

  return {

    fileSelected: (file)=>{

      var action = { type: 'FILE_SELECTED', data: file };
      dispatch(action);
    }
  };
}

const ListHeadContainer = ReactRedux.connect(

  mapStateToProps,
  mapDispatchToProps,
)(ListHead);

export default ListHeadContainer;