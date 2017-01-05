import * as ReactRedux from 'react-redux';
import ListBody from './list-body';

const mapStateToProps = (state, ownProps)=>{

  // state == redux store
  return {

    'userInfo': state['userInfo'],
    'fileData': state['fileData']
  };
}

const mapDispatchToProps = (dispatch, ownProps)=>{

  return {

    authorizeFile: (file)=>{

      var action = { type: 'AUTHORIZE_FILE', data: file };
      dispatch(action);
    }
  };
}

const ListBodyContainer = ReactRedux.connect(

  mapStateToProps,
  mapDispatchToProps,
)(ListBody);

export default ListBodyContainer;