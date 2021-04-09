import {defaultApiHelpers} from "./api";
import {Dispatch, MutableRefObject, useCallback, useEffect, useRef, useState} from "react";
import {VulcanState} from "../reducers/vulcan_reducer";
import {setSession, setStatus, VulcanAction} from "../actions/vulcan";
import {SessionStatusResponse, VulcanSession} from "../api_types";
import {Cancellable} from "etna-js/utils/cancellable";
import {hasNoRunningSteps} from '../selectors/workflow_selectors';

export const defaultSessionSyncHelpers = {
  statusIsFresh: true,
  requestPoll(post?: boolean) {
  },
};

function updateFromSessionResponse(response: SessionStatusResponse, state: MutableRefObject<VulcanState>, dispatch: Dispatch<VulcanAction>) {
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
  const [lastPollingRequest, setPollingRequest] = useState([false] as [boolean]);
  const [lastCompletedPollingRequest, setCompletedRequest] = useState(null as null | VulcanSession['inputs']);
  

  // We want to guarantee a minimum number of status polls,
  //   to make sure we get accurate status from the server.
  // Sometimes a status or submit POST will receive a response
  //   from the server where all steps show as PENDING,
  //   because the request just happened to catch the scheduler
  //   between tasks.
  // In that scenario, the UI will stop polling because there
  //   are no RUNNING steps, when in reality it just missed
  //   seeing a status update with RUNNING steps. By ensuring
  //   a set of status requests spread out over time, we
  //   catch this corner case and make the UI appear more responsive.
  const minimumNumberRequests = 3;
  const [currentRequestNumber, setCurrentRequestNumber] = useState(0 as number);

  // Note -- this will 'cancel' processing the result of previous requests when a new polling request is made, but
  // it does not actually cancel the underlying request.  Posts will thus complete just fine, although the result may
  // be superceeded by a requests poll.
  useEffect(() => {
    const [doPost] = lastPollingRequest;
    const cancellable = new Cancellable();

    if (!state.current.session.workflow_name) return;

    const sessionOfWork = state.current.session;
    const baseWork = doPost ? postInputs(sessionOfWork) : pollStatus(sessionOfWork);
    cancellable.race(baseWork).then(({result, cancelled}) => {
      if (cancelled || !result) return;
      updateFromSessionResponse(result, state, dispatch);
      setCompletedRequest(state.current.session.inputs);
    })

    // Reset our counter for each postInputs call
    if (doPost) setCurrentRequestNumber(0);

    return () => cancellable.cancel();
  }, [lastPollingRequest]);

  const requestPoll = useCallback((post = false) => {
    // Increment our request counter once all steps stop running
    if (hasNoRunningSteps(state.current.status)) {
      setCurrentRequestNumber(currentRequestNumber + 1);
    }

    if (!post &&
        hasNoRunningSteps(state.current.status) &&
        lastCompletedPollingRequest === state.current.session.inputs &&
        currentRequestNumber > minimumNumberRequests) {
      setCurrentRequestNumber(0);
      return;
    }

    setPollingRequest([post]);
  }, [lastCompletedPollingRequest]);

  useEffect(() => {
    const timerId = setInterval(() => requestPoll(), 1000);
    return () => clearInterval(timerId);
  }, [requestPoll])

  const statusIsFresh = lastCompletedPollingRequest === state.current.session.inputs;

  return {
    statusIsFresh,
    requestPoll,
  };
}