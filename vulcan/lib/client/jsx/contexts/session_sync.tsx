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
  const [lastPollingRequest, setPollingRequest] = useState([false, state.current.session.inputs] as [boolean, {[k: string]: any}]);

  // Note -- this will 'cancel' processing the result of previous requests when a new polling request is made, but
  // it does not actually cancel the underlying request.  Posts will thus complete just fine, although the result may
  // be superceeded by a requests poll.
  useEffect(() => {
    const [doPost] = lastPollingRequest;
    const cancellable = new Cancellable();

    if (!state.current.session.workflow_name) return;

    const baseWork = doPost ? postInputs(state.current.session) : pollStatus(state.current.session);
    cancellable.race(baseWork).then(({result, cancelled}) => {
      if (cancelled || !result) return;
      updateFromSessionResponse(result, state, dispatch);
    })

    return () => cancellable.cancel();
  }, [lastPollingRequest]);

  const requestPoll = useCallback((post = false) => {
    const [lastWasPost, lastInputs] = lastPollingRequest;

    if (lastInputs === state.current.session.inputs && lastWasPost) {
      return;
    }

    setPollingRequest([post, state.current.session.inputs]);
  }, []);

  useEffect(() => {
    const timerId = setInterval(() => requestPoll(), 1000);
    return () => clearInterval(timerId);
  }, [])

  return {
    requestPoll,
  };
}