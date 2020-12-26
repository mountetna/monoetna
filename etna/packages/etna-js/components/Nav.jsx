import React, {useEffect} from 'react';
require('./Nav.css');

import Icon from './icon';

const ICONS = {
  superuser: 'user-ninja',
  administrator: 'user-astronaut',
  editor: 'user-edit',
  viewer: 'user'
};

const Login = ({user, project_name}) => {
  if (!user) return null;

  let {first, last, permissions} = user;

  let role = (permissions[project_name] || {}).role;

  if (permissions.administration && permissions.administration.role == 'administrator') role = 'superuser';

  return (
    <div className='etna-login'>
      {first} {last}
      <Icon className='etna-user' icon={ ICONS[role] } title={role}/>
    </div>
  );
};

const Logo = ({LogoImage}) =>
  <div className='etna-logo'>
    <a href='/'>
      <LogoImage/>
    </a>
  </div>;

const Nav = ({logo, project_name, children, user}) => 
  <div className='etna-nav'>
    <Logo LogoImage={logo}/>
    {children}
    <Login user={user} project_name={project_name}/>
  </div>;

export default Nav;
