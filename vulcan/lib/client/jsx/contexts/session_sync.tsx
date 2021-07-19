import {defaultApiHelpers} from "./api";
import {Dispatch, MutableRefObject} from "react";
import {VulcanState} from "../reducers/vulcan_reducer";
import {finishPolling, setSession, setStatus, startPolling, VulcanAction} from "../actions/vulcan_actions";
import {SessionStatusResponse} from "../api_types";
import {hasNoRunningSteps} from "../selectors/workflow_selectors";
import {runPromise, useAsyncCallback} from "etna-js/utils/cancellable_helpers";

export const defaultSessionSyncHelpers = {
  requestPoll(post?: boolean): Promise<unknown> {
    return Promise.resolve();
  }, cancelPolling() {
  }
};

export function delay(ms: number) {
  return new Promise(resolve => setTimeout(resolve, ms));
}

function updateFromSessionResponse(response: SessionStatusResponse, dispatch: Dispatch<VulcanAction>) {
  dispatch(setStatus(response.status));
  dispatch(setSession(response.session));
}

export function useSessionSync(state: MutableRefObject<VulcanState>,
  showErrors: typeof defaultApiHelpers.showErrors,
  pollStatus: typeof defaultApiHelpers.pollStatus,
  postInputs: typeof defaultApiHelpers.postInputs,
  dispatch: Dispatch<VulcanAction>,
): typeof defaultSessionSyncHelpers {
  /*
     Begins the process of exchanging the current session with the server and updating our state based on this.
     For post requests, the inputs are sent and work is scheduled if it is not already.
     For non post requests, the inputs are sent but work is not scheduled.
     After that, it continues to poll until all steps are not running.
   */
  const [requestPoll, cancelPolling] = useAsyncCallback(function* (post = false) {
    dispatch(startPolling());
    if (!state.current.session.workflow_name) {
      return;
    }

    const baseWork = post ? postInputs(state.current.session) : pollStatus(state.current.session);
    const response = yield* runPromise(showErrors(baseWork));
    updateFromSessionResponse(response, dispatch);

    while (!hasNoRunningSteps(state.current.status)) {
      yield delay(3000);
      const response = yield* runPromise(showErrors(pollStatus(state.current.session)));
      updateFromSessionResponse(response, dispatch);
    }
  }, [dispatch, pollStatus, postInputs, state], () => {
    dispatch(finishPolling())
  });

  return {
    requestPoll, cancelPolling
  };
}