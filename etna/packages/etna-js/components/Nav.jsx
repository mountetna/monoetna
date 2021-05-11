import React, {useState} from 'react';
require('./Nav.css');
import { isSuperuser, isSuperEditor, isSuperViewer } from '../utils/janus';

import Icon from './icon';

const ICONS = {
  superuser: 'user-ninja',
  supereditor: 'user-nurse',
  superviewer: 'user-tie',
  administrator: 'user-astronaut',
  editor: 'user-edit',
  viewer: 'user'
};

const Login = ({user}) => {
  if (!user) return null;

  let {name, permissions} = user;

  let role = (permissions[CONFIG.project_name] || {}).role;

  if (isSuperuser(user)) role = 'superuser';
  else if (isSuperEditor(user)) role = 'supereditor';
  else if (isSuperViewer(user)) role = 'superviewer';

  return (
    <div className='etna-login' title={role}>
      <Icon className='etna-user' icon={ ICONS[role] }/>
      {name}
    </div>
  );
};

const Logo = ({LogoImage}) =>
  <div className='etna-logo'>
    <a href='/'>
      <LogoImage/>
    </a>
  </div>;

const Link = ({app}) => {
  let image = <img title={app} className='etna-link' src={ `/images/${app}.svg` }/>;

  let host_key = `${app}_host`;

  if (!CONFIG[host_key]) return image;

  let link = new URL(...[CONFIG.project_name, CONFIG[host_key]].filter(_=>_));

  return <a href={link}>{image}</a>;
}

const Control = ({currentApp}) => {
  let [ shown, setShown ] = useState(false);
  let apps = [ 'timur', 'metis', 'janus' ].filter(a => a != currentApp);

  if (!shown) return <div className='etna-control'>
    <div className='etna-control-show' onClick={ () => setShown(true) } >
      <Icon icon='bars'/>
    </div>
  </div>;

  return <div className='etna-control'>
    <div className='etna-control-links'>
      { apps.map( app => <Link key={app} app={app}/>) }
    </div>
    <div className='etna-control-hide' onClick={ () => setShown(false) } >
      <Icon icon='bars'/>
    </div>
  </div>;
}

const Nav = ({logo, app, children, user}) => 
  <nav className='etna-nav'>
    <Logo LogoImage={logo}/>
    {findValidChildren(children)}
    <Login user={user}/>
    <Control currentApp={app}/>
  </nav>

function findValidChildren(children) {
  if (children) {
    // Only return non-null children from an Array, or
    //  the singular child if provided. When only
    //  one child is passed in, it is an Object, not an Array,
    //  and has no `filter` method.
    if (Array.isArray(children)) {
      if (children.filter(_=>_).length) {
        return children;
      }
    } else {
      return children;
    }
  }
  return <div style={{flex: 1}}/>;
}

export default Nav;
