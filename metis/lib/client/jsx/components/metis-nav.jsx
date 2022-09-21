import * as React from 'react';
import Nav from 'etna-js/components/Nav';
import {selectUser} from 'etna-js/selectors/user-selector';
import {useReduxState} from 'etna-js/hooks/useReduxState';

const Logo = () => <div id='logo'/>;

const MetisNav = () => {
  let user = useReduxState( state => selectUser(state) );
  return <Nav logo={Logo} user={user} app='metis'/>;
};

export default MetisNav;
