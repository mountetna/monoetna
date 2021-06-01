import React, {useContext, useState, useEffect, useCallback, useMemo} from 'react';

import {VulcanContext} from '../../../contexts/vulcan_context';

import StepName from './step_name';

import {
  allDataNonEmpty, allExpectedOutputSources,
  dataOfSource,
  statusOfStep, stepOfSource,
} from "../../../selectors/workflow_selectors";
import {patchInputs} from "../../../actions/vulcan";
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
  }, [dispatch]);

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

  // Keeping the component `key` as just the step.name
  //   (without allStepsComplete)
  //   makes sure the components don't unmount / remount
  //   when they move from the Completed list to Pending.
  // This could happen if a user goes back to change
  //   a previous input step.
  // Allowing the component to unmount / remount means that
  //   the input values could reset to initial state and
  //   wipe out other elements of the user selection.
  // For example, in a multi-multiselect-string-all,
  //   removing one of three choices could force the other
  //   two to also reset if the component re-mounts.
  // This does also mean that the steps may not automatically
  //   hide their contents (i.e. the drawer close effect)
  //   if edited and Run is clicked, but that seems more minor
  //   compared to resetting user input values.
  const inner = 'isGroup' in step ?
      <StepUserInputDrawer
          key={`${step.name}`}
          step={step}
          handleInputChange={handleInputChange}
      /> :
      <StepUserInput
          key={`${step.name}`}
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
