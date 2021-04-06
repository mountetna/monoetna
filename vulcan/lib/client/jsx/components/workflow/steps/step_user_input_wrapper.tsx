import React, {useContext, useState, useEffect, useCallback, useMemo} from 'react';

import {VulcanContext} from '../../../contexts/vulcan_context';

import StepName from './step_name';

import {
  allDataNonEmpty, allExpectedOutputSources,
  dataOfSource,
  statusOfStep, stepOfSource,
} from "../../../selectors/workflow_selectors";
import {patchInputs, setInputs} from "../../../actions/vulcan";
import {UIStep} from "../user_interactions/inputs/input_types";
import StepUserInputDrawer from "./step_user_input_drawer";
import StepUserInput from "./step_user_input";
import {STATUS} from "../../../api_types";

export default function StepUserInputWrapper({step}: { step: UIStep['step'] }) {
  const [open, setOpen] = useState(true);
  const {dispatch, state} = useContext(VulcanContext);

  const {status, workflow, data, session} = state;

  const handleInputChange = useCallback((inputName: string, value: any) => {
    dispatch(patchInputs({[inputName]: value}));
  }, []);

  const toggleInputs = useCallback(() => setOpen(!open), [setOpen, open]);

  // Memoize this so that we can close or open the group when its effective value changes, rather than
  // also closing it every time any incidental change to data or session happens that doesn't result in
  // an effective change to this condition (see the useEffect below).
  const allInnerSources = useMemo(() => allExpectedOutputSources(step), [step]);
  const allInnerStepNames: string[] = useMemo(() => allInnerSources.map(stepOfSource).filter(v => v != null),
      [allInnerSources]) as string[];
  const allInnerStatus = useMemo(() => allInnerStepNames.map(stepName => statusOfStep(stepName, status)), [allInnerStepNames, status])
  const allStepsComplete = useMemo(() => allInnerStatus.every(s => s && s.status === STATUS.COMPLETE), [allInnerStatus]);
  const allDatas = useMemo(() => allInnerSources.map(source =>
      dataOfSource(source, workflow, status, data, session)), [allInnerSources, workflow, status, data, session]);
  const allDataIsNonEmpty = useMemo(() => allDataNonEmpty(allDatas), [allDatas]);

  useEffect(() => {
    if (!allDataIsNonEmpty) {
      setOpen(true);
    } else if (allStepsComplete) {
      setOpen(false);
    }
  }, [allDataIsNonEmpty, allStepsComplete]);

  const inner = 'isGroup' in step ?
      <StepUserInputDrawer
          key={`${step.name}-${allStepsComplete}`}
          step={step}
          handleInputChange={handleInputChange}
      /> :
      <StepUserInput
          key={`${step.name}-${allStepsComplete}`}
          step={step}
          handleInputChange={handleInputChange}
      />;

  return (
      <div className='step-user-input step'>
        <div onClick={toggleInputs}>
          <StepName step={step}/>
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
