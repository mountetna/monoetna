import * as ReactRedux from 'react-redux';
import ProjectEdit from './project-edit';

const mapStateToProps = (state, ownProps)=>{

  // state == redux store
  return {

    adminInfo: state['adminInfo']
  };
}

const mapDispatchToProps = (dispatch, ownProps)=>{

  return {};
}

const ProjectEditContainer = ReactRedux.connect(

  mapStateToProps,
  mapDispatchToProps,
)(ProjectEdit);

export default ProjectEditContainer;