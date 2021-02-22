import React, {useState, useContext, useEffect} from 'react';

import Icon from 'etna-js/components/icon';
import {VulcanContext} from '../../../contexts/vulcan';
import {STATUS} from '../../../models/steps';
import StepName from './step_name';

export default function Step({step, index}) {
  const {status, pathIndex, setStepIndex} = useContext(VulcanContext);

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
    <StepName
      step={step}
      status={stepStatus}
      onClick={handleOnClick}
    ></StepName>
  );
}
