import {Dispatch} from 'react';
import {setDownloadedData, VulcanAction} from '../actions/vulcan_actions';
import {defaultApiHelpers} from './api';
import {VulcanState} from '../reducers/vulcan_reducer';
import {shouldDownloadStep} from '../selectors/workflow_selectors';
import {useAsync} from 'etna-js/utils/cancellable_helpers';
import {runAttempts} from 'etna-js/utils/retryable';
import { MultiFileContentResponse, OutputConfig } from '../api_types';

export function useDataBuffering(
    state: VulcanState,
    dispatch: Dispatch<VulcanAction>,
    showErrors: typeof defaultApiHelpers.showErrors,
    readFiles: typeof defaultApiHelpers.readFiles,
) {
  const {status, workflow, data, projectName, workspaceId} = state;

  useAsync(function* () {
    if (!workspaceId) return;

    // Find the next status with a download that is an input to a ui-output or ui-query,
    // initiate a download, and let it run.
    for (let [step, stepStatus] of Object.entries(status.steps)) {
      if (stepStatus.status=="complete" && stepStatus.outputs?.files) {
        for (let outputName of stepStatus.outputs.files) {
          if (outputName in Object.keys(data)) continue;
          if (!shouldDownloadStep(stepStatus.name, workflow, outputName)) continue;

          const downloaded: MultiFileContentResponse = yield* runAttempts(() => readFiles(projectName, workspaceId, [outputName]));

          dispatch(setDownloadedData(outputName, downloaded[0]['content']));
          return;
        }
      }

    }
  }, [status, workflow, data]);
}
