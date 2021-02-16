import React, {useState, useContext, useEffect} from 'react';

import Icon from 'etna-js/components/icon';
import {VulcanContext} from '../../../contexts/vulcan';

// Hacked in for dev
const STATUS = [];
STATUS[0] = require('../../../../../server/data/status1.json');
STATUS[1] = require('../../../../../server/data/status2.json');
STATUS[2] = require('../../../../../server/data/status3.json');
STATUS[3] = require('../../../../../server/data/status4.json');
STATUS[4] = require('../../../../../server/data/status5.json');

// No status for step 6 yet...
STATUS[5] = require('../../../../../server/data/status5.json');

export default function Step({step, index, active}) {
  const {setStepIndex, setStatus} = useContext(VulcanContext);

  function handleOnClick(index) {
    setStepIndex(index);
    setStatus(STATUS[index]);
  }

  const icons = {
    complete: 'check',
    pending: 'clock',
    error: 'times-circle'
  };

  let icon = icons[step.status ? step.status : 'pending'];
  console.log(icon);
  return (
    <li className={active ? 'active step' : 'step'}>
      <div className='step-status-icon-wrapper'>
        <div className='step-status-icon-positioner'>
          <Icon className='step-status-icon' icon={icon}></Icon>
        </div>
      </div>
      {active ? (
        <div className='active-wrapper'>
          <div className='active-indicator'></div>
        </div>
      ) : null}
      <button onClick={() => handleOnClick(index)}>{step.name}</button>
    </li>
  );
}
