import * as React from 'react';

import TitleBar  from './nav/title-bar';
import MenuBarContainer   from './nav/menu-bar-container';

import UserEditContainer from './user-admin/user-edit-container';
import ProjectEditContainer from './user-admin/project-edit-container';
import PermissionEditContainer from './user-admin/permission-edit-container';

export default class UserAdminView extends React.Component{

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

          <img src='/img/logo_basic.png' alt='' />
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