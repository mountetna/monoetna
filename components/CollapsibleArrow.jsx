import React, { useCallback, useState } from 'react';
require('./CollapsibleArrow.css');

const onChangeDefault = () => null;
export default function CollapsibleArrow({ collapsed: collapsedProp, onChange: onChangeProp = onChangeDefault, children, label }) {
  const [collapsedState, setCollapsed] = useState(false);
  const collapsed = collapsedProp == null ? collapsedState : collapsedProp;

  const onChange = useCallback(() => {
    if (onChangeProp(!collapsed) !== false) {
      setCollapsed(!collapsed);
    }
  }, [onChangeProp, collapsed, setCollapsed]);

  if (collapsed) {
    return <div>
      <div className='etna-collapsible-arrow collapsed' onClick={onChange}/> {label}
    </div>;
  }

  return <div>
    <div>
      <div className='etna-collapsible-arrow' onClick={onChange}/> {label}
    </div>
    {children}
  </div>;
}
