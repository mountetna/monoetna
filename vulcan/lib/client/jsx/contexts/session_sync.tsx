import * as _ from 'lodash';

import {defaultApiHelpers} from "./api";
import {Dispatch, MutableRefObject, useCallback, useEffect, useRef, useState} from "react";
import {VulcanState} from "../reducers/vulcan_reducer";
import {setSession, setStatus, VulcanAction} from "../actions/vulcan";
import {SessionStatusResponse, VulcanSession} from "../api_types";
import {Cancellable} from "etna-js/utils/cancellable";
import {hasNoRunningSteps} from "../selectors/workflow_selectors";

export const defaultSessionSyncHelpers = {
  isPolling: false,
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
  const [curPolling, setCurPolling] = useState(0);
  const [_, setCancellable] = useState(new Cancellable());

  const requestPoll = useCallback((post = false) => {
    if (!post && hasNoRunningSteps(state.current.status)) {
      return;
    }

    if (!state.current.session.workflow_name) return;

    const cancellable = new Cancellable();
    setCancellable(c => {
      c.cancel();
      return cancellable;
    })

    const sessionOfWork = state.current.session;
    const baseWork = post ? postInputs(sessionOfWork) : pollStatus(sessionOfWork);

    function sync(work: Promise<SessionStatusResponse>): Promise<void> {
      return cancellable.race(work).then(({result, cancelled}) => {
        if (cancelled || !result) return;
        updateFromSessionResponse(result, dispatch);

        if (hasNoRunningSteps(state.current.status)) {
          return Promise.resolve();
        }

        return new Promise((resolve) => setTimeout(resolve, 1000))
          .then(() => sync(pollStatus(state.current.session)));
      });
    }

    setCurPolling(v => v + 1);
    sync(baseWork).then(r => console.log('finished polling')).finally(() => setCurPolling(v => v - 1));
  }, [dispatch, pollStatus, postInputs, state]);

  // Kind, of silly.  Just want this useEffect to run once on first mount.  But the linter will complain unless we
  // wrap the callback in a ref.  we do not want to fire this for every single change to the callback.
  const requestPollRef = useRef(requestPoll);
  useEffect(() => requestPollRef.current(), []);

  return {
    isPolling: curPolling > 0,
    requestPoll,
  };
}