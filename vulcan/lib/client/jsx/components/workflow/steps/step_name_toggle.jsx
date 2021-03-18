import React from 'react';

import StepName from './step_name';

export default function StepNameToggle({step, status}) {
  return (
    <div className='step-name-toggle-wrapper'>
      <div className='step-name-toggle arrow-toggle toggle-up-down'></div>
      <StepName step={step} status={status}></StepName>
    </div>
  );
}
