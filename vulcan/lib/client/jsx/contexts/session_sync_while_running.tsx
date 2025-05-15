import { defaultApiHelpers } from './api';
import {Dispatch, MutableRefObject} from 'react';
import {VulcanState} from '../reducers/vulcan_reducer';
import {checkCommittedStepPending, finishPolling, setRunId, useUIAccounting, setStatusFromStatuses, startPolling, VulcanAction, setLastConfig} from '../actions/vulcan_actions';
import {AccountingReturn, RunReturn, RunStatus, WorkspaceStatus} from '../api_types';
import { hasRunningSteps, paramValuesToRaw } from '../selectors/workflow_selectors';
import {runPromise, useAsyncCallback} from 'etna-js/utils/cancellable_helpers';
import {Maybe} from '../selectors/maybe';

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
  showErrors: typeof defaultApiHelpers.showErrors,
  requestRun: typeof defaultApiHelpers.requestRun,
  pullRunStatus: typeof defaultApiHelpers.pullRunStatus,
  getIsRunning: typeof defaultApiHelpers.getIsRunning,
  dispatch: Dispatch<VulcanAction>,
): typeof defaultSessionSyncHelpers {
  /*
    Initiates a hand-over of inputs' source of truth to the calculation server until either cancellation or there are no running jobs.

    If post == true: Ensures ui values were synced from session to server.
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
    const {projectName, workspaceId, status} = state;
    let {configId, runId} = state;
    let isRunning = false;

    dispatch(startPolling());
    if (!workspaceId) {
      console.log("Skipping post or work, no workspace_id")
      return;
    }

    if (startWork) {
      if (!configId) {
        console.log("Skipping work, no configId")
        return;
      }
      const response2: RunReturn = yield* runPromise(showErrors(requestRun(projectName,workspaceId,configId)));
      updateFromRunRequest(response2, dispatch);
      runId = response2.run_id;
      isRunning = true;
    } else {
      const runningResponse = yield* runPromise(showErrors(getIsRunning(projectName,workspaceId)));
      isRunning = runningResponse['running'];
    }

    while (isRunning) {
      yield delay(3000);
      if (!runId) {
        console.log("Skipping polling, no runId, ending polling")
        isRunning = true;
        return;
      }
      const runningResponse = yield* runPromise(showErrors(getIsRunning(projectName,workspaceId)));
      isRunning = runningResponse['running'];
      const response: RunStatus = yield* runPromise(showErrors(pullRunStatus(projectName,workspaceId,runId)));
      // Notably, after the first submission, changes in hash that occur do not have authority
      // to clear stale inputs.
      updateFromRunStatus(response, dispatch, isRunning);
    }
  }, [dispatch, requestRun, getIsRunning, pullRunStatus], () => {
    dispatch(finishPolling());
  });
  

  return {
    requestPoll, cancelPolling
  };
}
