import React, {useState, useContext, useEffect} from 'react';
import {VulcanContext} from '../../../contexts/vulcan';

export default function Step({step, index, onClick}) {
  const {stepIndex} = useContext(VulcanContext);

  return (
    <li className={stepIndex === index ? 'active-step' : ''}>{step.name}</li>
  );
}
