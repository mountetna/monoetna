import React, {useContext, useState, useEffect, useCallback, useMemo} from 'react';

import {VulcanContext} from '../../../contexts/vulcan_context';

import {STATUS} from '../../../api_types';
import UserInput from '../user_interactions/inputs/user_input';
import StepName from './step_name';

import {WorkflowStep} from "../../../api_types";
import {
    sourceNamesOfStep,
    statusOfStep,
    uiQueryOfStep, stepInputDataRaw,
} from "../../../selectors/workflow_selectors";
import {setInputs} from "../../../actions/vulcan";
import {InputSpecification, InputType, UIStep} from "../user_interactions/inputs/input_types";
import StepUserInputDrawer from "./step_user_input_drawer";
import StepUserInput from "./step_user_input";

export default function StepUserInputWrapper({step}: {step: UIStep['step']}) {
    const [open, setOpen] = useState(true);
    const {dispatch, state} = useContext(VulcanContext);

    const status = statusOfStep(step, state.status);
    const uiQuery = uiQueryOfStep(step);

    if (!status || !uiQuery) return null;

    const handleInputChange = useCallback((inputName: string, value: any) => {
        dispatch(setInputs({ ...state.inputs, [inputName]: value }));
    }, [state, dispatch, setInputs]);

    const toggleInputs = useCallback(() => setOpen(!open), [setOpen, open]);

    useEffect(() => {
        let stepStatus = status.status;
        if (STATUS.COMPLETE === stepStatus) {
            setOpen(false);
        } else if (STATUS.PENDING === stepStatus) {
            setOpen(true);
        }
    }, [status.status]);

    const inner = 'isGroup' in step ?
        <StepUserInputDrawer
            key={`${step.name}-${status.status}`}
            step={step}
            handleInputChange={handleInputChange}
        /> :
        <StepUserInput
            key={`${step.name}-${status.status}`}
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
