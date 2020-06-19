import * as React from 'react';

export default ({icon, className, overlay, title}) =>
  <span className='icon-grp fa-stack fa-fw' title={title}>
    <i className={ `icon ${className} ${icon} fa fa-2x fa-${icon}` }/>
    { overlay && <i className={ `overlay ${className} ${overlay} fa fa-stack-1x fa-${overlay}` }/> }
  </span>;
