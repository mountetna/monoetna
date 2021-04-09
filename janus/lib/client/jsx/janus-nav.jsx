import * as React from 'react';
import Nav from 'etna-js/components/Nav';
import {selectUser} from 'etna-js/selectors/user-selector';
import {isSuperViewer} from 'etna-js/utils/janus';
import {useReduxState} from 'etna-js/hooks/useReduxState';

const Logo = () => <div id='logo'/>;

const NavBar = ({user}) => <div id='nav'>
  <div className='nav_item'><a href='/settings'>Settings</a></div>
  { isSuperViewer(user) && <div className='nav_item'><a href='/admin'>Admin</a></div> }
</div>

const JanusNav = () => {
  let user = useReduxState( state => selectUser(state) );
  return <Nav logo={Logo} user={user} app='janus'>
    <NavBar user={user}/>
  </Nav>;
}

export default JanusNav;
