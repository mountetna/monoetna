import {defaultApiHelpers} from "./api";
import {Dispatch, MutableRefObject, useEffect, useRef} from "react";
import {VulcanState} from "../reducers/vulcan_reducer";
import {setSession, setStatus, VulcanAction} from "../actions/vulcan";
import * as _ from 'lodash';
import {SessionStatusResponse} from "../api_types";

export const defaultSessionSyncHelpers = {
    requestPoll() {
    },
};

export function useSessionSync(
    state: MutableRefObject<VulcanState>,
    scheduleWork: typeof defaultApiHelpers.scheduleWork,
    pollStatus: typeof defaultApiHelpers.pollStatus,
    postInputs: typeof defaultApiHelpers.postInputs,
    dispatch: Dispatch<VulcanAction>,
): typeof defaultSessionSyncHelpers {
    // 'switch' like mutex that indicates a desire for queueing a future poll/post quickly, while ensuring
    // only one is running at a given time.
    const loadingState = useRef({running: false, rid: 0, nextPollWithPost: false});

    useEffect(() => {
        const timerId = setInterval(() => requestPoll(), 3000);
        return () => clearInterval(timerId);
    }, [])

    function updateFromSessionResponse(response: SessionStatusResponse) {
        // Don't bother changing the state and trigger a bunch of other downstream stuff unless
        // the net effect is an actual change.
        if (!_.isEqual(state.current.status, response.status))
            dispatch(setStatus(response.status));
        if (!_.isEqual(state.current.session, response.session))
            dispatch(setSession(response.session));
    }

    function requestPoll(post = false) {
        const {running, rid} = loadingState.current;

        if (!running) {
            loadingState.current.nextPollWithPost = false;
            loadingState.current.running = true;
            const baseWork = post ? postInputs(state.current.session) : pollStatus(state.current.session);;
            scheduleWork(baseWork.then(updateFromSessionResponse)).finally(() => {
                loadingState.current.running = false;
                if (rid !== loadingState.current.rid) {
                    requestPoll(loadingState.current.nextPollWithPost);
                }
            });
        } else {
            loadingState.current.nextPollWithPost ||= post;
            loadingState.current.rid++;
        }
    }

    useEffect(() => {
        requestPoll();
    }, [state.current.session.inputs, state.current.workflow]);

    return {
        requestPoll,
    };
}