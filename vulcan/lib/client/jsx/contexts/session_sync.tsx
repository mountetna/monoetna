import {defaultApiHelpers} from "./api";
import {Dispatch, MutableRefObject, useCallback, useEffect, useRef, useState} from "react";
import {VulcanState} from "../reducers/vulcan_reducer";
import {finishPolling, setSession, setStatus, startPolling, VulcanAction} from "../actions/vulcan_actions";
import {SessionStatusResponse, VulcanSession} from "../api_types";
import {Cancellable} from "etna-js/utils/cancellable";
import {hasNoRunningSteps} from "../selectors/workflow_selectors";
import {showMessages} from "etna-js/actions/message_actions";
import {useActionInvoker} from "etna-js/hooks/useActionInvoker";

export const defaultSessionSyncHelpers = {
  requestPoll(post?: boolean) {
    return Promise.resolve();
  },
  cancelPolling() {},
  run() { return Promise.resolve(true); },
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
    invoke: ReturnType<typeof useActionInvoker>,
    validationErrors: VulcanState['validationErrors'],
    dispatch: Dispatch<VulcanAction>,
): typeof defaultSessionSyncHelpers {
  const [_, setCancellable] = useState(new Cancellable());

  const requestPoll = useCallback((post = false) => {
    // Guarantee that any call to requestPoll immediately cancels any existing one.
    const cancellable = new Cancellable();
    setCancellable(c => {
      c.cancel();
      return cancellable;
    })

    if (!post && hasNoRunningSteps(state.current.status)) {
      console.log('request to poll ignored, no running steps and not a post');
      return Promise.resolve();
    }

    if (!state.current.session.workflow_name) {
      console.log('current session does not have workflow_name set');
      return Promise.resolve();
    }

    const sessionOfWork = state.current.session;
    const baseWork = post ? postInputs(sessionOfWork) : pollStatus(sessionOfWork);

    function sync(work: Promise<SessionStatusResponse>, attempt=1): Promise<void> {
      return cancellable.race(work).then(({result, cancelled}) => {
        if (cancelled || !result) return;
        updateFromSessionResponse(result, dispatch);

        if (hasNoRunningSteps(state.current.status)) {
          console.log('running steps have completed, stopping polling');
          return Promise.resolve();
        }

        return new Promise((resolve) => setTimeout(resolve, 1000))
          .then(() => sync(pollStatus(state.current.session)));
      }, (e) => {
        // On failure, retry after a longer period of time.
        return new Promise((resolve) => setTimeout(resolve, Math.min(attempt + 1, 6) ** 2 * 1000))
          .then(() => sync(pollStatus(state.current.session), attempt + 1));
      });
    }

    dispatch(startPolling());
    return sync(baseWork).finally(() => dispatch(finishPolling()));
  }, [dispatch, pollStatus, postInputs, state]);

  const cancelPolling = useCallback(() => {
    setCancellable(c => {
      c.cancel();
      return new Cancellable();
    })
  }, [setCancellable]);

  const run = useCallback(() => {
    if (Object.keys(validationErrors).length > 0) {
      invoke(
        showMessages(
          Object.entries(validationErrors)
            .map(([inputName, validation]: [string, any]) => {
              let {
                inputLabel,
                errors
              }: {inputLabel: string; errors: string[]} = validation;
              return errors.map((e: string) => `${inputLabel}: ${e}`);
            })
            .flat()
        )
      );

      return Promise.resolve(false);
    } else {
      return requestPoll(true).then(() => true);
    }
  }, [requestPoll, invoke, validationErrors]);

  return {
    requestPoll,
    cancelPolling,
    run,
  };
}