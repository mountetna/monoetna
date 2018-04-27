import * as React from 'react';
import * as ReactRedux from 'react-redux';

import TitleBar  from './nav/title-bar';
import MenuBar   from './nav/menu-bar';
import ListHead  from './list/list-head';
import ListBody  from './list/list-body';
import LoginPanel from './auth/login-panel';

class MetisUI extends React.Component {
  constructor(){
    super();
  }

  renderLoginView(){
    return (
      <div id='listing-group'>
        <LoginPanel />
      </div>
    );
  }

  renderContent() {
    let { downloads, uploads, fails } = this.props.files;

    return (
      <div id='listing-group'>
        <table id='listing-table'>
          <ListHead />
          {(Object.keys(downloads).length || Object.keys(uploads).length || fails.length)?
            <ListBody /> : <tbody></tbody>
          }
        </table>
      </div>
    );
  }

  render(){
    return (
      <div id='metis-group'>
        <div id='header-group'>
          <TitleBar />
          <MenuBar />
        </div>
        <div className='logo-group'>
          <img src='/img/metis_logo_simple.png' alt='' />
        </div>
        <div id='left-column-group'>
        </div>
        { this.renderContent() }
      </div>
    );
  }
}

const MetisUIContainer = ReactRedux.connect(
  // map state
  ({user, files}) => ({user, files})
)(MetisUI);

export default MetisUIContainer;
