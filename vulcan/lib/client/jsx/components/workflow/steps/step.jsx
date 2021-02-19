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
    <li className={active ? 'active step' : 'step'}>
      <div className='step-status-icon-wrapper'>
        <div className='step-status-icon-positioner'>
          <Icon className={className} icon={icon.icon}></Icon>
        </div>
      </div>
      {active ? (
        <div className='active-wrapper'>
          <div className='active-indicator'></div>
        </div>
      ) : null}
      <button onClick={() => handleOnClick(index)}>
        {step.label || step.name}
      </button>
    </li>
  );
}
