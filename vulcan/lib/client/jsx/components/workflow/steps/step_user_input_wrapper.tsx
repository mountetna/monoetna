import React, {useContext, useState, useEffect, useCallback} from 'react';

import {VulcanContext} from '../../../contexts/vulcan_context';

import {STATUS} from '../../../api_types';
import StepName from './step_name';

import {
    allDataNonEmpty, allExpectedOutputSources,
    dataOfSource,
    statusOfStep,
    uiQueryOfStep,
} from "../../../selectors/workflow_selectors";
import {setInputs} from "../../../actions/vulcan";
import {UIStep} from "../user_interactions/inputs/input_types";
import StepUserInputDrawer from "./step_user_input_drawer";
import StepUserInput from "./step_user_input";

export default function StepUserInputWrapper({step}: { step: UIStep['step'] }) {
    const [open, setOpen] = useState(true);
    const {dispatch, state} = useContext(VulcanContext);

    const {status, workflow, data, session} = state;

    const stepStatus = statusOfStep(step, status);
    const uiQuery = uiQueryOfStep(step);

    if (!stepStatus || !uiQuery) return null;

    const handleInputChange = useCallback((inputName: string, value: any) => {
        dispatch(setInputs({...state.inputs, [inputName]: value}));
    }, [state, dispatch, setInputs]);

    const toggleInputs = useCallback(() => setOpen(!open), [setOpen, open]);
    const statusStr = stepStatus.status;

    useEffect(() => {
        // A group's "output" is the definition of values for all of its composite inputs.
        if (allDataNonEmpty(allExpectedOutputSources(step).map(source =>
            dataOfSource(source, workflow, status, data, session)))) {
            setOpen(true);
        } else if (STATUS.COMPLETE === statusStr) {
            setOpen(false);
        }
    }, [status, workflow, data, session, statusStr]);

    const inner = 'isGroup' in step ?
        <StepUserInputDrawer
            key={`${step.name}-${statusStr}`}
            step={step}
            handleInputChange={handleInputChange}
        /> :
        <StepUserInput
            key={`${step.name}-${statusStr}`}
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
