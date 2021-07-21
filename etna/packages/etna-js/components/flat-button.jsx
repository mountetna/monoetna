import * as React from 'react';
import Icon from './icon';
require('./flat-button.scss');

const FlatButton = ({icon, className, disabled, title, label, onClick}) => (
  <div title={title} className={ `flat-button ${className} ${disabled ? 'disabled' : ''}` } onClick={onClick}>
    { label && <span className='flat-label'>{label}</span> }
    { icon && <Icon disabled={disabled} className='flat-icon' icon={icon} /> }
  </div>
)

export default FlatButton;

