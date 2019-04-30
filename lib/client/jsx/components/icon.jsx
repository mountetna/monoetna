import * as React from 'react';

export default ({icon, overlay, title}) =>
  <span className='fa-stack fa-fw' title={title}>
    <i className={ `icon ${icon} fa fa-2x fa-${icon}` }/>
    { overlay && <i className={ `overlay ${overlay} fa fa-stack-1x fa-${overlay}` }/> }
  </span>;
