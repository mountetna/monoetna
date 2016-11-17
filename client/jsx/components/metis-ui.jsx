import * as React from 'react'

import TitleBar  from './nav/title-bar';
import MenuBar   from './nav/menu-bar';
import ListHeadContainer  from './list/list-head-container';
import ListBody  from './list/list-body';
import ListEmpty from './list/list-empty';
import LoginPanelContainer from './auth/login-panel-container';

export default class MetisUI extends React.Component{

  constructor(){

    super();
  }

  listBody(fileList, fileUploads){

    return <ListBody fileList={ fileList } fileUploads={ fileUploads } />;
  }

  renderLoginView(){

    if(this['props']['appState']['loginStatus']){

      return <div>{ 'sup' }</div>;
    }
    else{

      return <LoginPanelContainer />;
    }
  }

  render(){

    var fileList = this['props']['appState']['fileList'];
    var fileUploads = this['props']['appState']['fileUploads'];

    return (

      <div id='metis-group'>
        
        <div id='header-group'>
          
          <TitleBar />
          <MenuBar />
        </div>
        <div className='logo-group'>

          <img src='/img/logo_dna_color_round.png' alt='' />
        </div>
        <div id='left-column-group'>
        </div>
        <div id='listing-group'>

          { this.renderLoginView() }
          {/*
          <ListHeadContainer />
          { (fileList.length || fileUploads.length) ? this.listBody(fileList, fileUploads) : <ListEmpty /> }
          */}
        </div>
      </div>
    );
  }
}