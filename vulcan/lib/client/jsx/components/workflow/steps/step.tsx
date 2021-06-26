import React from 'react';

import StepName from './step_name';
import {WorkflowStep} from "../../../api_types";

export default function Step({step}: {step: WorkflowStep}) {
  return (
    <div className='step'>
      <StepName step={step}/>
    </div>
  );
}
