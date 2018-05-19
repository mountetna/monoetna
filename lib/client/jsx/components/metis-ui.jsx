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

    let columnWidths = {
      type: '70px',
      name: '60%',
      updated: '30%',
      size: '10%',
      control: '100px'
    };

    return (
      <div id='listing-group'>
        <div id='listing-table'>
          <ListHead widths={ columnWidths } />
          {(Object.keys(downloads).length || Object.keys(uploads).length || fails.length)?
            <ListBody widths={ columnWidths }/> : <div/>
          }
        </div>
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
