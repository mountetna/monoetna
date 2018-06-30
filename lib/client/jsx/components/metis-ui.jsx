import * as React from 'react';
import { connect } from 'react-redux';

import TitleBar  from './nav/title-bar';
import MenuBar   from './nav/menu-bar';
import ListHead  from './list/list-head';
import ListBody  from './list/list-body';
import FolderBreadcrumb from './folder-breadcrumb';
import ControlBar from './control-bar';
import LoginPanel from './auth/login-panel';

class MetisUI extends React.Component {
  renderLoginView(){
    return (
      <div id='listing-group'>
        <LoginPanel />
      </div>
    );
  }

  componentDidMount() {
    let path = window.location.pathname;

    let [ _, folder_name ] = path.match(/browse\/(.*)$/) || [];

    this.props.setCurrentFolder(folder_name == undefined ? folder_name : decodeURI(folder_name));
    this.props.retrieveFiles(folder_name);
  }

  renderContent() {
    let columnWidths = {
      type: '90px',
      name: '60%',
      updated: '30%',
      size: '10%',
      control: '100px'
    };

    return (
      <div id='listing-group'>
        <ListHead widths={ columnWidths } />
        <ListBody widths={ columnWidths }/>
      </div>
    );
  }

  render() {
    return (
      <div id='metis-group'>
        <div id='header-group'>
          <TitleBar />
          <MenuBar />
        </div>
        <div id='logo-group'>
          <img src='/img/metis_logo_simple.png' alt='' />
        </div>
        <div id='control-group'>
          <FolderBreadcrumb/>
          <ControlBar/>
        </div>
        { this.renderContent() }
      </div>
    );
  }
}

const MetisUIContainer = connect(
  // map state
  null,
  // map dispatch
  (dispatch) => ({
    retrieveFiles: (folder_name) => dispatch({type: 'RETRIEVE_FILES', folder_name}),
    setCurrentFolder: (folder_name) => dispatch({type: 'SET_CURRENT_FOLDER', folder_name})
  })
)(MetisUI);

export default MetisUIContainer;
