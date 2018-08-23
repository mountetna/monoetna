import * as React from 'react';
import { connect } from 'react-redux';

const MenuBar = ({user}) =>
  <div id='nav-menu'>
    <div className='nav-user'>
      { user.first } { user.last }
    </div>
  </div>;

export default connect(
  // map state
  ({user}) => ({user})
)(MenuBar);
