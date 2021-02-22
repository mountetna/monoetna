import React, {useState, useContext, useEffect} from 'react';

import {VulcanContext} from '../../../contexts/vulcan';

import {validStep} from '../../../selectors/workflow_selector';

export default function StepPending({step, stepIndex}) {
  return <div className='step-pending'>Waiting for the step to finish.</div>;
}
