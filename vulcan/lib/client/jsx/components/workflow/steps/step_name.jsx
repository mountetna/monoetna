import React from 'react';

import Icon from 'etna-js/components/icon';
import {STATUS} from '../../../models/steps';

export default function StepName({step, status, onClick}) {
  const icons = {};
  icons[STATUS.COMPLETE] = {
    icon: 'check',
    className: 'light green'
  };
  icons[STATUS.PENDING] = {icon: 'clock', className: 'light'};
  icons[STATUS.ERROR] = {icon: 'times-circle', className: 'light red'};

  let icon = icons[status || STATUS.PENDING];

  let className = `step-status-icon ${icon.className}`;

  return (
    <div className='step-name'>
      <div className='step-status-icon-wrapper'>
        <Icon className={className} icon={icon.icon}></Icon>
      </div>
      <div className='step-button' onClick={onClick}>
        {step.label || step.name}
      </div>
    </div>
  );
}
