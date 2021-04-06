import {defaultApiHelpers} from "./api";
import {Dispatch, MutableRefObject, useCallback, useEffect, useRef, useState} from "react";
import {VulcanState} from "../reducers/vulcan_reducer";
import {setSession, setStatus, VulcanAction} from "../actions/vulcan";
import * as _ from 'lodash';
import {SessionStatusResponse, VulcanSession} from "../api_types";
import {Cancellable} from "etna-js/utils/cancellable";

export const defaultSessionSyncHelpers = {
  statusIsFresh: true,
  requestPoll(post?: boolean) {
  },
};

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
  const [lastCompletedPollingRequest, setCompletedRequest] = useState(null as null | VulcanSession['inputs']);

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
      setCompletedRequest(result.session.inputs);
    })

    return () => cancellable.cancel();
  }, [lastPollingRequest]);

  const requestPoll = useCallback((post = false) => {
    if (!post && state.current.status[0].every(s => s.status !== "running") && lastCompletedPollingRequest === state.current.session.inputs) {
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