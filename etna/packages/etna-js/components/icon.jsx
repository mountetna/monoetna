import * as React from 'react';

export default ({icon, className, disabled, overlay, title, onClick}) => (
  <span className={`icon-grp ${className ? `icon-grp-${className}` : ''} ${disabled ? 'disabled' : ''} fa-stack fa-fw`} title={title} onClick={disabled ? null : onClick}>
    <i className={`icon ${className} ${disabled ? 'disabled' : ''} ${icon} fa fa-2x fa-${icon}`} />
    {overlay && (
      <i
        className={`overlay ${className} ${disabled ? 'disabled' : ''} ${overlay} fa fa-stack-1x fa-${overlay}`}
      />
    )}
  </span>
);
