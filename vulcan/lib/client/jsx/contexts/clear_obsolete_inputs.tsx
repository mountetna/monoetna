import {VulcanState} from "../reducers/vulcan_reducer";
import {Dispatch, useEffect, useRef} from "react";
import * as _ from 'lodash';
import {setInputs, VulcanAction} from "../actions/vulcan";
import {SessionStatusResponse} from "../api_types";

export function useClearObsoleteInputs(state: VulcanState, dispatch: Dispatch<VulcanAction>) {
    const {inputs, workflow} = state;
    const lastInputsRef = useRef({inputs});

    useEffect(() => {
        if (workflow == null) return;
        const {dependencies_of_outputs} = workflow;

        const {current: {inputs: lastInputs}} = lastInputsRef;

        let newInputs: SessionStatusResponse['session']['inputs'] = inputs;

        function clearDependentInputs(source: string) {
            const dependencies = dependencies_of_outputs[source];
            if (dependencies) {
                dependencies.forEach(source => {
                    if (source in newInputs) {
                        if (newInputs === inputs) newInputs = {...inputs};
                        delete newInputs[source];
                    }
                });
            }
        }

        Object.keys(inputs).forEach(source => {
            if (source in lastInputs && !_.isEqual(inputs[source], lastInputs[source])) {
                clearDependentInputs(source);
            }
        });

        if (newInputs !== inputs) {
            dispatch(setInputs(newInputs));
        }

        lastInputsRef.current.inputs = inputs;
    }, [inputs, workflow]);
}