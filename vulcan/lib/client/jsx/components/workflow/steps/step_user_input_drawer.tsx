import React, {useContext, useMemo} from 'react';

import {VulcanContext} from '../../../contexts/vulcan_context';
import UserInput from '../user_interactions/inputs/user_input';
import {GroupedInputStep, InputSpecification} from "../user_interactions/inputs/input_types";
import {stepInputDataRaw, stepOfSource, stepOfStatus, uiQueryOfStep} from "../../../selectors/workflow_selectors";

export default function StepUserInputDrawer({step, handleInputChange}: {step: GroupedInputStep, handleInputChange: (sourceName: string, val: any) => void}) {
    let {state} = useContext(VulcanContext);
    const {workflow, session, status, data, inputs} = state;

    const uiQuery = uiQueryOfStep(step);

    if (!uiQuery) return null;
    if (!workflow) return null;

    // We need to unpack the grouped steps and add docs
    let stepInputs: InputSpecification[] = useMemo(() => step.in.map((input) => {
            // We need to fetch the original step, to see if the options
            //   data is available.
            const originalStepName = stepOfSource(input.source);
            const originalStep = originalStepName ? stepOfStatus(originalStepName, workflow) : null;
            const inputData = originalStep ? stepInputDataRaw(originalStep, status, data, session) : {};

            return {
                type: uiQuery,
                label: input.label || input.source,
                default: inputs[input.source] || null,
                data: inputData,
                name: input.source,
                doc: input.doc,
            };
        }),
[]);

    return (
        <React.Fragment>
            {stepInputs.map((input, index) => {
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