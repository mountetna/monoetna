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

  renderContent(){

    if(!this['props']['userInfo']['loginStatus']){

      return this.renderLoginView();
    }
    else{

      return this.renderRegular();
    }
  }

  renderLoginView(){

    return (

      <div id='listing-group'>
        
        <LoginPanelContainer />
      </div>
    );
  }

  renderRegular(){

    var fileList = this['props']['fileData']['fileList'];
    var fileUploads = this['props']['fileData']['fileUploads'];

    return (

      <div id='listing-group'>
        
        <table id='listing-table'>
        
          <ListHeadContainer />
          { (fileList.length || fileUploads.length) ? <ListBodyContainer /> : '' }
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

          <img src='/img/logo_dna_color_round.png' alt='' />
        </div>
        <div id='left-column-group'>
        </div>
        { this.renderContent() }
      </div>
    );
  }
}