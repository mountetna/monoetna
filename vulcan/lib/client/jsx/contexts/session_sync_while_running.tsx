import { defaultApiHelpers } from './api';
import {Dispatch, MutableRefObject} from 'react';
import {VulcanState} from '../reducers/vulcan_reducer';
import {checkCommittedStepPending, finishPolling, setRunId, useUIAccounting, setStatusFromStatuses, startPolling, VulcanAction} from '../actions/vulcan_actions';
import {AccountingReturn, RunReturn, RunStatus} from '../api_types';
import { hasRunningSteps } from '../selectors/workflow_selectors';
import {runPromise, useAsyncCallback} from 'etna-js/utils/cancellable_helpers';
import {Maybe} from '../selectors/maybe';

export const defaultSessionSyncHelpers = {
  requestPoll(post = false, submittingStep: Maybe<string> = null): Promise<unknown> {
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
  submittingStep: Maybe<string>
) {
  dispatch(useUIAccounting(response, submittingStep));
  dispatch(checkCommittedStepPending(submittingStep));
}

function updateFromRunRequest(
  response: RunReturn,
  dispatch: Dispatch<VulcanAction>,
  submittingStep: Maybe<string>,
) {
  dispatch(setRunId(response.run_id));
}

function updateFromRunStatus(
  response: RunStatus,
  dispatch: Dispatch<VulcanAction>,
  submittingStep: Maybe<string>,
) {
  dispatch(setStatusFromStatuses(response));
  dispatch(checkCommittedStepPending(submittingStep));
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
    Initiates a hand-over of inputs' source of truth to the calculation server until either cancellation or there are no running job.
    If post == true: Ensures ui values were synced from session to server, then begins the running of work.
    Afterwards or if post != true: requests step status updates regularly from server until all work completes.
    Note: only step statuses are updated automatically.  New file outputs are not retrieved here.
    
    Old: Begins the process of exchanging the current session with the server and updating our state based on this.
    For post requests, the inputs are sent *and work is scheduled* if it is not already.
    For non post requests, the inputs are sent but work is not scheduled.
    After that, it continues to poll until all steps are not running.
  */
  const {projectName, workspace, workspaceId, configId, status} = state.current;
  let {runId} = state.current;
  
  const [requestPoll, cancelPolling] = useAsyncCallback(function* (
    post = false,
    submittingStep: Maybe<string> = null,
  ) {
    dispatch(startPolling());
    if (!workspace || !workspaceId || !configId) {
      return;
    }

    // initial work
    if (post) {
      const response1: AccountingReturn = yield* runPromise(showErrors(postUIValues(projectName,workspaceId,status,submittingStep)));
      updateFromPostUIAccounting(response1, dispatch, submittingStep)
      // ToDo: Error if configId changed?
      const response2: RunReturn = yield* runPromise(showErrors(requestRun(projectName,workspaceId,response1.config_id)));
      updateFromRunRequest(response2, dispatch, submittingStep);
      runId = response2.run_id;
    }

    while (hasRunningSteps(state.current.status)) {
      yield delay(3000);
      if (!runId) return;
      const response: RunStatus = yield* runPromise(showErrors(pullRunStatus(projectName,workspaceId,runId)));
      // Notably, after the first submission, changes in hash that occur do not have authority
      // to clear stale inputs.
      updateFromRunStatus(response, dispatch, null);
    }
  }, [dispatch, pullRunStatus, postUIValues, state], () => {
    dispatch(finishPolling());
  });
  

  return {
    requestPoll, cancelPolling
  };
}
