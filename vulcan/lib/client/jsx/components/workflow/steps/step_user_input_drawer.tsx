import React, {useContext, useMemo} from 'react';

import {VulcanContext} from '../../../contexts/vulcan_context';
import UserInput from '../user_interactions/inputs/user_input';
import {
  GroupedInputStep,
  InputSpecification
} from '../user_interactions/inputs/input_types';
import {
  stepInputDataRaw,
  stepOfSource,
  stepOfStatus,
  uiQueryOfStep,
  sortInputsByLabel
} from '../../../selectors/workflow_selectors';

export default function StepUserInputDrawer({
  step,
  handleInputChange
}: {
  step: GroupedInputStep;
  handleInputChange: (sourceName: string, val: any) => void;
}) {
  let {state} = useContext(VulcanContext);
  const {workflow, session, status, data, inputs} = state;

  // We need to unpack the grouped steps and add docs
  let stepInputs: InputSpecification[] = useMemo(() => !workflow ? [] : step.in.map((input) => {
        // We need to fetch the original step, to see if the options are available and to unpack the true uiQuery
        const originalStepName = stepOfSource(input.source);
        const originalStep = originalStepName ? stepOfStatus(originalStepName, workflow) : null;
        const inputData = originalStep ? stepInputDataRaw(originalStep, status, data, session) : {};
        const uiQuery = originalStep ? uiQueryOfStep(originalStep) : null;

            return {
              type: uiQuery || 'default',
              label: input.label || input.source,
              value: inputs[input.source] || null,
              data: inputData,
              name: input.source,
              doc: input.doc
            };
          }),
    [step, workflow, status, data, session, inputs]
  );

  if (!workflow) return null;

  if (!workflow) return null;

  return (
    <React.Fragment>
      {sortInputsByLabel(stepInputs).map((input, index) => {
        return (
          <UserInput
            input={input}
            hideLabel={false}
            onChange={handleInputChange}
            key={index}
          />
        );
      })}
    </React.Fragment>
  );
}
