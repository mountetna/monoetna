import {defaultApiHelpers} from "./api";
import {Dispatch, MutableRefObject, useCallback, useEffect, useRef, useState} from "react";
import {VulcanState} from "../reducers/vulcan_reducer";
import {setSession, setStatus, VulcanAction} from "../actions/vulcan";
import * as _ from 'lodash';
import {SessionStatusResponse} from "../api_types";
import {Cancellable} from "etna-js/utils/cancellable";

export const defaultSessionSyncHelpers = {
    requestPoll(post?: boolean) {
    },
};

// Note -- this uses
function updateFromSessionResponse(response: SessionStatusResponse, state: MutableRefObject<VulcanState>, dispatch: Dispatch<VulcanAction>) {
    // Don't bother changing the state and trigger a bunch of other downstream stuff unless
    // the net effect is an actual change.
    if (!_.isEqual(state.current.status, response.status))
        dispatch(setStatus(response.status));
    if (!_.isEqual(state.current.session, response.session))
        dispatch(setSession(response.session));
}


export function useSessionSync(
    state: MutableRefObject<VulcanState>,
    scheduleWork: typeof defaultApiHelpers.scheduleWork,
    pollStatus: typeof defaultApiHelpers.pollStatus,
    postInputs: typeof defaultApiHelpers.postInputs,
    dispatch: Dispatch<VulcanAction>,
): typeof defaultSessionSyncHelpers {
    const [lastPollingRequest, setPollingRequest] = useState([false] as [boolean]);

    useEffect(() => {
        const timerId = setInterval(() => requestPoll(), 3000);
        return () => clearInterval(timerId);
    }, [])

    useEffect(() => {
        const [doPost] = lastPollingRequest;
        const cancellable = new Cancellable();

        const baseWork = doPost ? postInputs(state.current.session) : pollStatus(state.current.session);
        cancellable.race(baseWork).then(({result, cancelled}) => {
            if (cancelled || !result) return;
            return updateFromSessionResponse(result, state, dispatch);
        })

        return () => cancellable.cancel();
    }, [lastPollingRequest]);

    const requestPoll = useCallback((post = false) => {
      setPollingRequest([post]);
    }, []);

    return {
        requestPoll,
    };
}