import React, {useReducer, createContext, useEffect, useRef} from 'react';
import VulcanReducer, {defaultVulcanState, VulcanState} from '../reducers/vulcan_reducer';
import {VulcanAction} from "../actions/vulcan";
import {defaultSessionStorageHelpers, useLocalSessionStorage} from "./session_storage";
import {defaultApiHelpers, useApi} from "./api";
import {useWorkflowsLoading} from "./workflows_loading";
import {useDataBuffering} from "./data_buffering";

export const defaultContext = {
    state: defaultVulcanState as VulcanState,
    stateRef: {current: defaultVulcanState},
    // This would be set with a context dispatch when the provide is actually installed.
    dispatch: (a: VulcanAction) => console.log('action dispatched but not handled', a),
    ...defaultSessionStorageHelpers,
    ...defaultApiHelpers,
}

export type VulcanContextData = typeof defaultContext;
export const VulcanContext = createContext(defaultContext);
export type VulcanContext = typeof VulcanContext;

export const VulcanProvider = (props: {params: {}, state?: VulcanState, storage?: typeof localStorage}) => {
    const [state, dispatch] = useReducer(VulcanReducer, {...defaultVulcanState, ...(props.state || {})});
    const stateRef = useRef(state);
    const {getLocalSession} = useLocalSessionStorage(stateRef, props);
    const {isLoading, scheduleWork} = useApi(stateRef, props);
    useDataBuffering(state, dispatch, scheduleWork);
    useWorkflowsLoading(JSON.stringify(props.params), dispatch, scheduleWork);

    return (
        <VulcanContext.Provider value={{
            state,
            stateRef,
            dispatch,
            getLocalSession,
            isLoading,
            scheduleWork
        }}>
        </VulcanContext.Provider>
    );
};
