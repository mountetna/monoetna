import * as React from 'react';

import TitleBar  from './nav/title-bar';
import MenuBarContainer   from './nav/menu-bar-container';

import UserEditContainer from './admin/user-edit-container';
import ProjectEditContainer from './admin/project-edit-container';
import PermissionEditContainer from './admin/permission-edit-container';

export default class AdminUI extends React.Component{

  constructor(){

    super();
  }

  render(){

    var userInfo = this['props']['userInfo'];

    return (

      <div id='admin-group'>

        <div id='header-group'>
          
          <TitleBar />
          <MenuBarContainer />
        </div>
        <div className='logo-group'>

          <img src='/img/logo_dna_color_round.png' alt='' />
        </div>
        <div id='left-column-group'>
        </div>
        { (userInfo['masterPerms']) ? <UserEditContainer /> : '' }
        { (userInfo['masterPerms']) ? <ProjectEditContainer /> : '' }
        { (userInfo['masterPerms']) ? <PermissionEditContainer /> : '' }
      </div>
    );
  }
}