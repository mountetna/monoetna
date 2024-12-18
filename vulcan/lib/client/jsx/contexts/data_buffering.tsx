import {Dispatch} from 'react';
import {setFileContent, VulcanAction} from '../actions/vulcan_actions';
import {defaultApiHelpers} from './api';
import {VulcanState} from '../reducers/vulcan_reducer';
import {parseIfCan, shouldDownloadOutput} from '../selectors/workflow_selectors';
import {useAsync} from 'etna-js/utils/cancellable_helpers';
import {runAttempts} from 'etna-js/utils/retryable';
import { MultiFileContentResponse } from '../api_types';

export function useDataBuffering(
    state: VulcanState,
    dispatch: Dispatch<VulcanAction>,
    showErrors: typeof defaultApiHelpers.showErrors,
    readFiles: typeof defaultApiHelpers.readFiles,
) {
  const {status, workspace, projectName, workspaceId} = state;

  useAsync(function* () {
    if (!workspaceId || !workspace) return;

    // Find the next status with a download that is an input to a ui-output or ui-query,
    // initiate a download, and let it run.
    for (let [step, stepStatus] of Object.entries(status.steps)) {
      // ToDo: update this next line around where full outputs knowledge is kept!
      if (stepStatus.status=="complete" && !!workspace.steps[step].output?.files) {
        for (let outputName of workspace.steps[step].output.files) {
          if (outputName in Object.keys(status.file_contents)) continue;
          if (!shouldDownloadOutput(workspace, outputName)) continue;

          const downloaded: MultiFileContentResponse = yield* runAttempts(() => readFiles(projectName, workspaceId, [outputName]));

          dispatch(setFileContent(outputName, parseIfCan(downloaded[0]['content'])));
          return;
        }
      }

    }
  }, [status, workspace, projectName, workspaceId]);
}
