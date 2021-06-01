import {defaultApiHelpers} from "./api";
import {Dispatch, MutableRefObject, useCallback, useEffect, useRef, useState} from "react";
import {VulcanState} from "../reducers/vulcan_reducer";
import {setSession, setStatus, VulcanAction} from "../actions/vulcan";
import {SessionStatusResponse, VulcanSession} from "../api_types";
import {Cancellable} from "etna-js/utils/cancellable";
import {hasNoRunningSteps} from "../selectors/workflow_selectors";

export const defaultSessionSyncHelpers = {
  statusIsFresh: true,
  requestPoll(post?: boolean) {
  },
};

function updateFromSessionResponse(response: SessionStatusResponse, dispatch: Dispatch<VulcanAction>) {
  dispatch(setStatus(response.status));
  dispatch(setSession(response.session));
}

export function useSessionSync(
    state: MutableRefObject<VulcanState>,
    scheduleWork: typeof defaultApiHelpers.scheduleWork,
    pollStatus: typeof defaultApiHelpers.pollStatus,
    postInputs: typeof defaultApiHelpers.postInputs,
    dispatch: Dispatch<VulcanAction>,
): typeof defaultSessionSyncHelpers {
  const [lastPollingRequest, setPollingRequest] = useState([false, state.current] as [boolean, VulcanState]);
  const [lastCompletedPollingRequest, setCompletedRequest] = useState(null as null | VulcanSession['inputs']);

  // Note -- this will 'cancel' processing the result of previous requests when a new polling request is made, but
  // it does not actually cancel the underlying request.  Posts will thus complete just fine, although the result may
  // be superceeded by a requests poll.
  useEffect(() => {
    const [doPost, state] = lastPollingRequest;
    const cancellable = new Cancellable();

    if (!state.session.workflow_name) return;

    const sessionOfWork = state.session;
    const baseWork = doPost ? postInputs(sessionOfWork) : pollStatus(sessionOfWork);
    cancellable.race(baseWork).then(({result, cancelled}) => {
      if (cancelled || !result) return;
      updateFromSessionResponse(result, dispatch);
      setCompletedRequest(state.session.inputs);
    })

    return () => cancellable.cancel();
  }, [lastPollingRequest, dispatch, pollStatus, postInputs, state]);

  const statusIsFresh = lastCompletedPollingRequest === state.current.session.inputs;
  const requestPoll = useCallback((post = false) => {
    if (!post && (hasNoRunningSteps(state.current.status) || statusIsFresh)) {
      return;
    }

    setPollingRequest([post, { ...state.current }]);
  }, [state, statusIsFresh]);

  useEffect(() => {
    const timerId = setInterval(() => requestPoll(), 1000);
    return () => clearInterval(timerId);
  }, [requestPoll])

  // Kind, of silly.  Just want this useEffect to run once on first mount.  But the linter will complain unless we
  // wrap the callback in a ref.  we do not want to fire this for every single change to the callback.
  const requestPollRef = useRef(requestPoll);
  useEffect(() => requestPollRef.current(), []);

  return {
    statusIsFresh,
    requestPoll,
  };
}