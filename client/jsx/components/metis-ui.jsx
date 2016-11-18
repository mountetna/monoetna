import * as React from 'react'

import TitleBar  from './nav/title-bar';
import MenuBarContainer   from './nav/menu-bar-container';
import ListHeadContainer  from './list/list-head-container';
import ListBodyContainer  from './list/list-body-container';
import ListEmpty from './list/list-empty';
import LoginPanelContainer from './auth/login-panel-container';

export default class MetisUI extends React.Component{

  constructor(){

    super();
  }

  renderLoginView(){

    if(!this['props']['appState']['loginStatus']){

      return (

        <div id='listing-group'>
          
          <LoginPanelContainer />
        </div>
      );
    }
    else{

      var appState = this['props']['appState'];
      var fileList = appState['fileList'];
      var fileUploads = appState['fileUploads'];

      return (

        <div id='listing-group'>
          
          <ListHeadContainer />
          { (fileList.length || fileUploads.length) ? <ListBodyContainer /> : <ListEmpty /> }
        </div>
      );
    }
  }

  render(){
    
    return (

      <div id='metis-group'>
        
        <div id='header-group'>
          
          <TitleBar />
          <MenuBarContainer />
        </div>
        <div className='logo-group'>

          <img src='/img/logo_dna_color_round.png' alt='' />
        </div>
        <div id='left-column-group'>
        </div>
        { this.renderLoginView() }
      </div>
    );
  }
}