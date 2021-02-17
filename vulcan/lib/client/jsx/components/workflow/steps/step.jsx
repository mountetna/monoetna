import React, {useState, useContext, useEffect} from 'react';

import Icon from 'etna-js/components/icon';
import {VulcanContext} from '../../../contexts/vulcan';
import {STATUS} from '../../../models/steps';

// Hacked in for dev
const MOCK_STATUS = [];
MOCK_STATUS[0] = require('../../../../../server/data/status1.json');
MOCK_STATUS[1] = require('../../../../../server/data/status2.json');
MOCK_STATUS[2] = require('../../../../../server/data/status3.json');
MOCK_STATUS[3] = require('../../../../../server/data/status4.json');
MOCK_STATUS[4] = require('../../../../../server/data/status5.json');

// No status for step 6 yet...
MOCK_STATUS[5] = require('../../../../../server/data/status5.json');

export default function Step({step, index, active}) {
  const {setStepIndex, setStatus} = useContext(VulcanContext);

  function handleOnClick(index) {
    setStepIndex(index);
    setStatus(MOCK_STATUS[index]);
  }

  const icons = {};
  icons[STATUS.COMPLETE] = {
    icon: 'check',
    className: 'light green'
  };
  icons[STATUS.PENDING] = {icon: 'clock', className: 'light'};
  icons[STATUS.ERROR] = {icon: 'times-circle', className: 'light red'};

  let icon = icons[step.status || STATUS.PENDING];
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
