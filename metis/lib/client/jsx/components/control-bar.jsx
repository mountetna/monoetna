import * as React from 'react';
import { connect } from 'react-redux';
import Icon from 'etna-js/components/icon';
import { selectUserRole } from '../selectors/user-selector';

// this is the Metis control bar, which contains basic operations
// like 'upload file' and 'create folder'

const ControlButton = ({onClick, icon, overlay, title}) => {
  return <button className='control-btn' onClick={ onClick } title={ title }>
    <Icon icon={icon} overlay={overlay}/>
  </button>;
};

const ControlBar = ({children, buttons, user_role}) => {
  buttons = buttons.filter(({role: button_role}) => !button_role || user_role <= button_role);

  if (!buttons.length) return null;

  return <div id='control-bar'>
    { children }
    {
      buttons.map(button => <ControlButton key={ button.title } {...button}/>)
    }
  </div>;
};

export default connect(
  state => ({ user_role: selectUserRole(state) })
)(ControlBar);
