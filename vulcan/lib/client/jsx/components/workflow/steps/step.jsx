import React, {useContext} from 'react';

import {VulcanContext} from '../../../contexts/vulcan';
import StepName from './step_name';
import {STATUS} from '../../../models/steps';

export default function Step({step, index, onClick}) {
  const {status, pathIndex} = useContext(VulcanContext);

  if (null == pathIndex) return null;

  function handleOnClick() {
    onClick(index);
  }

  const stepStatus = status ? status[pathIndex][index].status : STATUS.PENDING;

  return (
    <div className='step'>
      <StepName
        step={step}
        status={stepStatus}
        onClick={handleOnClick}
      ></StepName>
    </div>
  );
}
