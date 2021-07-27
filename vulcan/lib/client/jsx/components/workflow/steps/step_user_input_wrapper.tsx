import React, {
  useContext, useState, useEffect, useCallback, useMemo
} from 'react';

import {VulcanContext} from '../../../contexts/vulcan_context';
import StepName from './step_name';
import { statusOfStep } from '../../../selectors/workflow_selectors';
import StepUserInputDrawer from './step_user_input_drawer';
import {STATUS} from '../../../api_types';
import {WorkflowStepGroup} from "../user_interactions/inputs/input_types";

export default function StepUserInputWrapper({group}: { group: WorkflowStepGroup }) {
  const [open, setOpen] = useState(true);
  const {state} = useContext(VulcanContext);

  const {status} = state;

  const toggleInputs = useCallback(() => setOpen(!open), [setOpen, open]);
  const allInnerStatus = useMemo(() => group.steps.map((step) => statusOfStep(step, status)), [group, status]);
  const allStepsComplete = useMemo(() => allInnerStatus.every((s) => s && s.status === STATUS.COMPLETE),
    [allInnerStatus]
  );

  const hasValidationErrors = group.steps.some((step) => state.validationErrors.some(([stepName]) => stepName === step.name));
  const shouldOpen = hasValidationErrors;
  const shouldClose = allStepsComplete && !hasValidationErrors;

  useEffect(() => {
    if (shouldClose) {
      setOpen(false);
    } else if (shouldOpen) {
      setOpen(true);
    }
  }, [shouldOpen, shouldClose]);

  return (<div
      className={`step-user-input step ${hasValidationErrors ? 'error' : ''}`}
    >
      <div onClick={toggleInputs}>
        <StepName step={group} showToggle={true} open={open}/>
      </div>
      <div
        className={`step-user-input-inputs sliding-panel vertical ${open ? 'open' : 'closed'}`}
      >
        <StepUserInputDrawer group={group}/>
      </div>
    </div>);
}
