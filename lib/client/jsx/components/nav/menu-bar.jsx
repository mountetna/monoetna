import * as React from 'react';
import { connect } from 'react-redux';
import { selectUserName, selectUserRole } from '../../selectors/user-selector';
import Icon from '../icon';

const ICONS = {
  administrator: 'user-astronaut',
  editor: 'user-edit',
  viewer: 'user'
};

const MenuBar = ({first, last, role, permissions}) => {
  console.log("role are");
  console.log(role);
  return <div id='nav-menu'>
    <div className='nav-user'>
      { first } { last }
      <Icon icon={ ICONS[role] } title={role}/>
    </div>
  </div>;
}

export default connect(
  // map state
  (state) => ({
    ...selectUserName(state),
    role: selectUserRole(state)
  })
)(MenuBar);
