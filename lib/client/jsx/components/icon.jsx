import * as React from 'react';

export default ({icon, overlay}) =>
  <span className='fa-stack fa-fw'>
    <i className={ `icon ${icon} fa fa-2x fa-${icon}` }/>
    { overlay && <i className={ `overlay ${overlay} fa fa-stack-1x fa-${overlay}` }/> }
  </span>;
