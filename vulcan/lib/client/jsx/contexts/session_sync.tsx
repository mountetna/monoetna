import {defaultApiHelpers} from "./api";
import {Dispatch, MutableRefObject, useCallback, useEffect, useRef, useState} from "react";
import {VulcanState} from "../reducers/vulcan_reducer";
import {finishPolling, setSession, setStatus, startPolling, VulcanAction} from "../actions/vulcan_actions";
import {SessionStatusResponse, VulcanSession} from "../api_types";
import {Cancellable} from "etna-js/utils/cancellable";
import {hasNoRunningSteps} from "../selectors/workflow_selectors";

export const defaultSessionSyncHelpers = {
  requestPoll(post?: boolean) {
  },
  cancelPolling() {}
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
  const [_, setCancellable] = useState(new Cancellable());

  const requestPoll = useCallback((post = false) => {
    if (!post && hasNoRunningSteps(state.current.status)) {
      console.log('request to poll ignored, no running steps and not a post');
      return;
    }

    if (!state.current.session.workflow_name) {
      console.log('current session does not have workflow_name set');
      return;
    }

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
          console.log('running steps have completed, stopping polling');
          return Promise.resolve();
        }

        return new Promise((resolve) => setTimeout(resolve, 1000))
          .then(() => sync(pollStatus(state.current.session)));
      });
    }

    dispatch(startPolling());
    sync(baseWork).finally(() => dispatch(finishPolling()));
  }, [dispatch, pollStatus, postInputs, state]);

  const cancelPolling = useCallback(() => {
    setCancellable(c => {
      c.cancel();
      return new Cancellable();
    })
  }, [setCancellable]);

  return {
    requestPoll,
    cancelPolling,
  };
}