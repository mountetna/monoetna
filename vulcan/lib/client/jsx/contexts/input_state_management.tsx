import {VulcanState} from '../reducers/vulcan_reducer';
import {Dispatch, useEffect} from 'react';
import {
  patchInputs,
  removeDownloads,
  removeInputs,
  VulcanAction
} from '../actions/vulcan_actions';
import {
  isPendingUiQuery,
  missingOutputsForStep,
  pendingSteps,
  sourceNameOfReference,
  statusOfStep,
  stepOfSource
} from '../selectors/workflow_selectors';

export function useInputStateManagement(
  state: VulcanState,
  dispatch: Dispatch<VulcanAction>,
  isPolling: boolean
) {
}
