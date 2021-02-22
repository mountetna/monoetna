import React, {useState, useContext, useEffect} from 'react';

import Icon from 'etna-js/components/icon';
import {VulcanContext} from '../../../contexts/vulcan';
import {STATUS} from '../../../models/steps';

export default function Step({step, index, active}) {
  const {status, pathIndex, setStepIndex, setStatus} = useContext(
    VulcanContext
  );

  function handleOnClick(index) {
    setStepIndex(index);
  }

  if (!status || !pathIndex) return null;

  const stepStatus = status[pathIndex][index].status;

  const icons = {};
  icons[STATUS.COMPLETE] = {
    icon: 'check',
    className: 'light green'
  };
  icons[STATUS.PENDING] = {icon: 'clock', className: 'light'};
  icons[STATUS.ERROR] = {icon: 'times-circle', className: 'light red'};

  let icon = icons[stepStatus || STATUS.PENDING];

  let className = `step-status-icon ${icon.className}`;

  return (
    <div className={active ? 'active step' : 'step'}>
      <div className='step-status-icon-wrapper'>
        <Icon className={className} icon={icon.icon}></Icon>
      </div>
      <div className='step-button' onClick={() => handleOnClick(index)}>
        {step.label || step.name}
      </div>
      <div className='step-active'>
        {active && <Icon className='active-indicator' icon='caret-left'/> }
      </div>
    </div>
  );
}
