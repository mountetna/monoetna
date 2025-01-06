import { defaultApiHelpers } from './api';
import {Dispatch, MutableRefObject} from 'react';
import {VulcanState} from '../reducers/vulcan_reducer';
import {checkCommittedStepPending, finishPolling, setRunId, useUIAccounting, setStatusFromStatuses, startPolling, VulcanAction, setLastConfig} from '../actions/vulcan_actions';
import {AccountingReturn, RunReturn, RunStatus, WorkspaceStatus} from '../api_types';
import { hasRunningSteps, paramValuesToRaw } from '../selectors/workflow_selectors';
import {runPromise, useAsyncCallback} from 'etna-js/utils/cancellable_helpers';
import {Maybe} from '../selectors/maybe';

export const defaultSessionSyncHelpers = {
  requestPoll(post = false, startWork = false, submittingStep?: string | null): Promise<unknown> {
    return Promise.resolve();
  }, cancelPolling() {
  }
};

export function delay(ms: number) {
  return new Promise(resolve => setTimeout(resolve, ms));
}

function updateFromPostUIAccounting(
  response: AccountingReturn,
  dispatch: Dispatch<VulcanAction>,
  status: WorkspaceStatus,
  submittingStep: string | null
) {
  if (submittingStep === null) {
    dispatch(setLastConfig(paramValuesToRaw(status.params)))
  }
  dispatch(useUIAccounting(response, submittingStep));
}

function updateFromRunRequest(
  response: RunReturn,
  dispatch: Dispatch<VulcanAction>,
) {
  dispatch(setRunId(response.run_id));
}

function updateFromRunStatus(
  response: RunStatus,
  dispatch: Dispatch<VulcanAction>
) {
  dispatch(setStatusFromStatuses(response));
}

export function useSessionSyncWhileRunning(
  state: MutableRefObject<VulcanState>,
  showErrors: typeof defaultApiHelpers.showErrors,
  requestRun: typeof defaultApiHelpers.requestRun,
  pullRunStatus: typeof defaultApiHelpers.pullRunStatus,
  postUIValues: typeof defaultApiHelpers.postUIValues,
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
  const {projectName, workspace, workspaceId, status} = state.current;
  let {configId, runId} = state.current;
  
  const [requestPoll, cancelPolling] = useAsyncCallback(function* (
    post = false,
    startWork = false,
    submittingStep?: string | null,
  ) {
    dispatch(startPolling());
    if (!workspace || !workspaceId || !configId) {
      return;
    }

    // First: Post UI selections
    if (post) {
      if (submittingStep===undefined) {
        console.error('Cannot request posting of inputs without a submittingStep, bug in client');
        return;
      }
      const response1: AccountingReturn = yield* runPromise(showErrors(postUIValues(projectName,workspaceId,status,submittingStep)));
      updateFromPostUIAccounting(response1, dispatch,status,submittingStep)
      configId = response1.config_id;
    }

    if (startWork) {
      const response2: RunReturn = yield* runPromise(showErrors(requestRun(projectName,workspaceId,configId)));
      updateFromRunRequest(response2, dispatch);
      runId = response2.run_id;
    }

    while (hasRunningSteps(state.current.status)) {
      yield delay(3000);
      if (!runId) return;
      const response: RunStatus = yield* runPromise(showErrors(pullRunStatus(projectName,workspaceId,runId)));
      // Notably, after the first submission, changes in hash that occur do not have authority
      // to clear stale inputs.
      updateFromRunStatus(response, dispatch);
    }
  }, [dispatch, pullRunStatus, postUIValues, state], () => {
    dispatch(finishPolling());
  });
  

  return {
    requestPoll, cancelPolling
  };
}
