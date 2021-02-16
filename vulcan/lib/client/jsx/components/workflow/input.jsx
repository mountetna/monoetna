import React, {useState, useContext, useEffect} from 'react';

import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
import {showMessages} from 'etna-js/actions/message_actions';

import {VulcanContext} from '../../contexts/vulcan';
import {getWorkflow, getWorkflows} from '../../api/vulcan';
import StepInput from './steps/step_input';

export default function Input() {
  const invoke = useActionInvoker();
  const {
    workflow,
    pathIndex,
    stepIndex,
    setPathIndex,
    setStepIndex
  } = useContext(VulcanContext);

  return (
    <section className='step-input'>
      <StepInput></StepInput>
    </section>
  );
}
