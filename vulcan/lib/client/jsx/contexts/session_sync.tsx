import {defaultApiHelpers} from "./api";
import {Dispatch, MutableRefObject, useCallback, useEffect, useRef, useState} from "react";
import {VulcanState} from "../reducers/vulcan_reducer";
import {finishPolling, setSession, setStatus, startPolling, VulcanAction} from "../actions/vulcan_actions";
import {SessionStatusResponse} from "../api_types";
import {Cancellable, CancelOrResult} from "etna-js/utils/cancellable";
import {hasNoRunningSteps} from "../selectors/workflow_selectors";

export const defaultSessionSyncHelpers = {
  requestPoll(post?: boolean): [Promise<CancelOrResult<SessionStatusResponse>>, Promise<boolean>] {
    return [Promise.resolve({cancelled: true}), Promise.resolve(false)];
  }, cancelPolling() {
  }
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

  /*
     Begins the process of exchanging the current session with the server and updating our state based on this.
     1.  Cancels any existing polling process, ensuring that no race condition allows for an existing polling process
         to update state after a call to requestPoll.
     2.  post=true => the server should begin to process any runnable steps, where as false will only report back the status.
     3.  Returns two promises -- the first, results from the first syncing of the session to the server, containing updated
         hashes based on those inputs.  The second is a continuous (retry enabled) process to poll the server so long as
         any status='running' process exists.
   */
  const requestPoll = useCallback<(post: boolean) => [Promise<CancelOrResult<SessionStatusResponse>>, Promise<boolean>]>(
    (post = false) => {
      // Guarantee that any call to requestPoll immediately cancels any existing one.
      const cancellable = new Cancellable();
      setCancellable(c => {
        c.cancel();
        return cancellable;
      })

      if (!state.current.session.workflow_name) {
        console.log('current session does not have workflow_name set');
        return [Promise.resolve({cancelled: true}), Promise.resolve(true)];
      }

      const sessionOfWork = state.current.session;
      const baseWork = post ? postInputs(sessionOfWork) : pollStatus(sessionOfWork);

      function sync(work: Promise<SessionStatusResponse>, attempt = 1): Promise<CancelOrResult<SessionStatusResponse>> {
        return cancellable.race(work).then(({result, cancelled}) => {
          if (!cancelled && result) updateFromSessionResponse(result, dispatch);
          return {result, cancelled};
        }, (e) => {
          // On failure, retry after a longer period of time.
          return new Promise((resolve) => setTimeout(resolve, Math.min(attempt + 1, 6) ** 2 * 1000))
            .then(() => sync(pollStatus(state.current.session), attempt + 1));
        });
      }

      function continuePolling({result, cancelled}: CancelOrResult<SessionStatusResponse>): Promise<boolean> {
        if (cancelled || !result) return Promise.resolve(false);

        if (hasNoRunningSteps(result.status)) {
          console.log('running steps have completed, stopping polling');
          return Promise.resolve(true);
        }

        return new Promise((resolve) => setTimeout(resolve, 1000))
          .then(() => sync(pollStatus(result.session)).then(continuePolling));
      }

      dispatch(startPolling());
      const firstSync = sync(baseWork);

      return [firstSync, firstSync.then(continuePolling).finally(() => dispatch(finishPolling()))];
    },
    [dispatch, pollStatus, postInputs, state]
  );

  const cancelPolling = useCallback(() => {
    setCancellable(c => {
      c.cancel();
      return new Cancellable();
    })
  }, [setCancellable]);

  return {
    requestPoll, cancelPolling
  };
}