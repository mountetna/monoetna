import { defaultApiHelpers } from './api';
import {Dispatch, MutableRefObject} from 'react';
import {VulcanState} from '../reducers/vulcan_reducer';
import {checkCommittedStepPending, finishPolling, setRunId, useUIAccounting, setStatusFromStatuses, startPolling, VulcanAction, setLastConfig} from '../actions/vulcan_actions';
import {AccountingReturn, RunReturn, RunStatus, WorkspaceStatus} from '../api_types';
import { hasRunningSteps, paramValuesToRaw } from '../selectors/workflow_selectors';
import {runPromise, useAsyncCallback} from 'etna-js/utils/cancellable_helpers';
import {Maybe} from '../selectors/maybe';
import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
import {showMessages} from 'etna-js/actions/message_actions';

export const defaultSessionSyncHelpers = {
  requestPoll(state: VulcanState, startWork = false): Promise<unknown> {
    return Promise.resolve();
  }, cancelPolling() {
  }
};

export function delay(ms: number) {
  return new Promise(resolve => setTimeout(resolve, ms));
}

function updateFromRunRequest(
  response: RunReturn,
  dispatch: Dispatch<VulcanAction>,
) {
  dispatch(setRunId(response.run_id));
}

function updateFromRunStatus(
  response: RunStatus,
  dispatch: Dispatch<VulcanAction>,
  isRunning: boolean,
) {
  dispatch(setStatusFromStatuses(response, isRunning));
}

export function useSessionSyncWhileRunning(
  invoke: ReturnType<typeof useActionInvoker>,
  showErrors: typeof defaultApiHelpers.showErrors,
  requestRun: typeof defaultApiHelpers.requestRun,
  pullRunStatus: typeof defaultApiHelpers.pullRunStatus,
  getIsRunning: typeof defaultApiHelpers.getIsRunning,
  dispatch: Dispatch<VulcanAction>,
): typeof defaultSessionSyncHelpers {
  /*
    Initiates a hand-over of workspaces' source of truth to the calculation server until either cancellation or there are no running jobs.

    Submission of inputs is handled elsewhere now.
    If startwork == true: Begins the scheduling of work.
    Afterwards / always: requests step status updates regularly from server until all work completes.
    Note: only step statuses are updated automatically.  New file outputs are not retrieved here.
    
    ---
    Old: Begins the process of exchanging the current session with the server and updating our state based on this.
    For post requests, the inputs are sent *and work is scheduled* if it is not already.
    For non post requests, the inputs are sent but work is not scheduled.
    After that, it continues to poll until all steps are not running.
  */
  
  const [requestPoll, cancelPolling] = useAsyncCallback(function* (
    state: VulcanState,
    startWork = false,
  ) {
    const {projectName, workspaceId, configId} = state;
    let {runId, isRunning} = state;
    const showError = (e: string) => {
      // invoke(dismissMessages());
      console.error(`likely a bug, please alert Dan! requestPoll Error: ${e}`)
      invoke(showMessages([`Error: ${e}`]))
    };

    if (!workspaceId) {
      showError("Skipping polling request, workspaceId unknown")
      return;
    }

    if (startWork) {
      if (isRunning) {
        showError("Skipping new 'Run' request and awaiting what is already running to finish");
      } else {
        if (!configId) {
          showError("Skipping 'Run' request, configId unknown");
          return;
        }
        const runResponse: RunReturn = yield* runPromise(showErrors(requestRun(projectName,workspaceId,configId)));
        updateFromRunRequest(runResponse, dispatch);
        runId = runResponse.run_id;
      }
    }

    isRunning = true;
    while (isRunning) {
      yield delay(3000);
      if (runId==null) {
        showError("Polling run status not possible, runId unknown")
        return;
      }
      const runningResponse = yield* runPromise(showErrors(getIsRunning(projectName,workspaceId)));
      isRunning = runningResponse['running'];
      const response: RunStatus = yield* runPromise(showErrors(pullRunStatus(projectName,workspaceId,runId)));
      updateFromRunStatus(response, dispatch, isRunning);
    }
  }, [dispatch, requestRun, getIsRunning, pullRunStatus], () => {
    // ToDo: cancelRunning!
  });

  return {
    requestPoll, cancelPolling
  };
}
