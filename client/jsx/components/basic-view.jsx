import * as React from 'react';

import TitleBar  from './nav/title-bar';
import MenuBarContainer   from './nav/menu-bar-container';
import LoginPanelContainer from './auth/login-panel-container';

export default class BasicView extends React.Component{

  constructor(){

    super();
  }

  renderContent(){

    if(!this['props']['userInfo']['loginStatus']){

      return this.renderLoginView();
    }
  }

  renderLoginView(){

    return (

      <div id='listing-group'>
        
        <LoginPanelContainer />
      </div>
    );
  }

  render(){

    return (

      <div id='polyphemus-group'>
        
        <div id='header-group'>
          
          <TitleBar />
          <MenuBarContainer />
        </div>
        <div className='logo-group'>

          { /*<img src='/img/logo_dna_color_round.png' alt='' />*/}
        </div>
        <div id='left-column-group'>
        </div>
        { this.renderContent() }
      </div>
    )
  }
}