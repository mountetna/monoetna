import * as React from 'react';
import { connect } from 'react-redux';
import { selectUserName, selectUserRole } from '../selectors/user-selector';
import Icon from './icon';

const ICONS = {
  administrator: 'user-astronaut',
  editor: 'user-edit',
  viewer: 'user'
};

const MetisNav = ({first, last, role, permissions}) =>
  <div id='metis-nav'>
    <div id='logo-group'>
      <div id='logo'/>
    </div>
    <div id='nav-user'>
      { first } { last }
      <Icon icon={ ICONS[role] } title={role}/>
    </div>
  </div>

export default connect(
  // map state
  (state) => ({
    ...selectUserName(state),
    role: selectUserRole(state)
  })
)(MetisNav);
