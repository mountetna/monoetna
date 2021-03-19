import React, {useContext} from 'react';

import {VulcanContext} from '../../../contexts/vulcan';
import StepName from './step_name';
import {STATUS} from '../../../models/steps';

export default function Step({step, index}) {
  const {status, pathIndex} = useContext(VulcanContext);

  if (null == pathIndex) return null;

  const stepStatus =
    status && status[pathIndex] && status[pathIndex][index]
      ? status[pathIndex][index].status
      : STATUS.PENDING;

  return (
    <div className='step'>
      <StepName step={step} status={stepStatus}></StepName>
    </div>
  );
}
