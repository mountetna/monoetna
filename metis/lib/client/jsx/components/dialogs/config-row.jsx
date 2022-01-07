import * as React from 'react';

const ConfigRow = ({children, label}) => (
  <div className='config-row'>
    <div className='label'>{label}</div>
    <div className='input'>{children}</div>
  </div>
);

export default ConfigRow;
