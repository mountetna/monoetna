import React, {useState, useContext, useEffect} from 'react';
import {VulcanContext} from '../../../contexts/vulcan';

export default function Step({step, index, active}) {
  const {setStepIndex} = useContext(VulcanContext);

  function handleOnClick(index) {
    setStepIndex(index);
  }

  return (
    <li className={active ? 'active step' : 'step'}>
      {active ? (
        <div className='active-wrapper'>
          <div className='active-indicator'></div>
        </div>
      ) : null}
      <button onClick={() => handleOnClick(index)}>{step.name}</button>
    </li>
  );
}
