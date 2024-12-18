import React from 'react';

import StepName from './step_name';
import {WorkspaceStep} from '../../../api_types';

export default function Step({step}: {step: WorkspaceStep}) {
  return (
    <div className='step'>
      <StepName step={step}/>
    </div>
  );
}
