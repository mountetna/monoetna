import * as React from 'react';

import TitleBar  from './nav/title-bar';
import MenuBarContainer   from './nav/menu-bar-container';
import ListHeadContainer  from './list/list-head-container';
import ListBodyContainer  from './list/list-body-container';
import LoginPanelContainer from './auth/login-panel-container';

export default class MetisUI extends React.Component{
  constructor(){
    super();
  }

  renderLoginView(){
    return (
      <div id='listing-group'>
        <LoginPanelContainer />
      </div>
    );
  }

  renderContent() {
    let { downloads, uploads, fails } = this.props.files;

    return (
      <div id='listing-group'>
        <table id='listing-table'>
          <ListHeadContainer />
          {(Object.keys(downloads).length || Object.keys(uploads).length || fails.length)?
            <ListBodyContainer /> : <tbody></tbody>
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
          <MenuBarContainer />
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
