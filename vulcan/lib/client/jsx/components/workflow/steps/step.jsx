import React, {useContext} from 'react';

import {VulcanContext} from '../../../contexts/vulcan';
import StepName from './step_name';

export default function Step({step, index, onClick}) {
  const {status, pathIndex} = useContext(VulcanContext);

<<<<<<< HEAD
  if (!status || null === pathIndex) return null;
=======
  if (!status || null == pathIndex) return null;
>>>>>>> cs/vulcan-error-handling

  function handleOnClick() {
    onClick(index);
  }

  const stepStatus = status[pathIndex][index].status;

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
