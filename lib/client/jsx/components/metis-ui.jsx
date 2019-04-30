import * as React from 'react';
import { connect } from 'react-redux';
import { selectCurrentFolder } from '../selectors/directory-selector';

import MetisNav from './metis-nav';
import FolderView from './folder-view';
import ModalDialog from './modal-dialog';

class MetisUI extends React.Component {
  componentDidMount() {
    let path = window.location.pathname;

    let [ _, folder_name ] = path.match(/browse\/(.*)$/) || [];

    this.props.setCurrentFolder(folder_name == undefined ? folder_name : decodeURI(folder_name));
    this.props.retrieveFiles(folder_name);
  }


  render() {
    let { current_folder } = this.props;
    return (
      <div id='metis-group'>
        <MetisNav/>
        <FolderView folder={current_folder}/>
        <ModalDialog/>
      </div>
    );
  }
}

const MetisUIContainer = connect(
  // map state
  (state) => ({current_folder: selectCurrentFolder(state)}),

  // map dispatch
  (dispatch) => ({
    retrieveFiles: (folder_name) => dispatch({type: 'RETRIEVE_FILES', folder_name}),
    setCurrentFolder: (folder_name) => dispatch({type: 'SET_CURRENT_FOLDER', folder_name})
  })
)(MetisUI);

export default MetisUIContainer;
