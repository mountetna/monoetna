import React, {
  useContext,
  useState,
  useEffect,
  useCallback,
  useMemo
} from 'react';

import {VulcanContext} from '../../../contexts/vulcan_context';

import StepName from './step_name';

import {
  allDataNonEmpty,
  allExpectedOutputSources,
  dataOfSource,
  statusOfStep,
  stepOfSource,
  sourceNamesOfStep
} from '../../../selectors/workflow_selectors';
import {patchInputs} from '../../../actions/vulcan';
import {UIStep} from '../user_interactions/inputs/input_types';
import StepUserInputDrawer from './step_user_input_drawer';
import StepUserInput from './step_user_input';
import {STATUS} from '../../../api_types';

export default function StepUserInputWrapper({step}: {step: UIStep['step']}) {
  const [open, setOpen] = useState(true);
  const {dispatch, state} = useContext(VulcanContext);

  const {status, workflow, data, session} = state;

  const handleInputChange = useCallback(
    (inputName: string, value: any) => {
      dispatch(patchInputs({[inputName]: value}));
    },
    [dispatch]
  );

  const toggleInputs = useCallback(() => setOpen(!open), [setOpen, open]);

  let inner;
  let inputNames: string[] = [];

  if ('isGroup' in step) {
    inner = (
      <StepUserInputDrawer
        key={`${step.name}`}
        step={step}
        handleInputChange={handleInputChange}
      />
    );
    inputNames = step.in.map((i) => i.source);
  } else {
    inner = (
      <StepUserInput
        key={`${step.name}`}
        step={step}
        handleInputChange={handleInputChange}
      />
    );
    inputNames = sourceNamesOfStep(step).map((outputName) => outputName);
  }

  // Memoize this so that we can close or open the group when its effective value changes, rather than
  // also closing it every time any incidental change to data or session happens that doesn't result in
  // an effective change to this condition (see the useEffect below).
  const allInnerSources = useMemo(() => allExpectedOutputSources(step), [step]);
  const allInnerStepNames: string[] = useMemo(
    () => allInnerSources.map(stepOfSource).filter((v) => v != null),
    [allInnerSources]
  ) as string[];
  const allInnerStatus = useMemo(
    () => allInnerStepNames.map((stepName) => statusOfStep(stepName, status)),
    [allInnerStepNames, status]
  );
  const allStepsComplete = useMemo(
    () => allInnerStatus.every((s) => s && s.status === STATUS.COMPLETE),
    [allInnerStatus]
  );
  const allDatas = useMemo(
    () =>
      allInnerSources.map((source) =>
        dataOfSource(source, workflow, status, data, session)
      ),
    [allInnerSources, workflow, status, data, session]
  );
  const allDataIsNonEmpty = useMemo(() => allDataNonEmpty(allDatas), [
    allDatas
  ]);
  const hasValidationErrors = inputNames.some((inputName) =>
    state.validationErrors.hasOwnProperty(inputName)
  );

  useEffect(() => {
    if (!allDataIsNonEmpty || hasValidationErrors) {
      setOpen(true);
    } else if (allStepsComplete && !hasValidationErrors) {
      setOpen(false);
    }
  }, [allDataIsNonEmpty, allStepsComplete, hasValidationErrors]);

  return (
    <div
      className={`step-user-input step ${hasValidationErrors ? 'error' : ''}`}
    >
      <div onClick={toggleInputs}>
        <StepName step={step} />
      </div>
      <div
        className={`step-user-input-inputs sliding-panel vertical ${
          open ? 'open' : 'closed'
        }`}
      >
        {inner}
      </div>
    </div>
  );
}
