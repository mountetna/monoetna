import React, {useState, useContext, useEffect} from 'react';

export default function StepError({step}) {
  return (
    <div className='step-error'>
      {step.error || 'Something went wrong with this step.'}
    </div>
  );
}
