import {Dispatch, useEffect} from 'react';
import {endSyncing, startSyncing, updateFiles, useUIAccounting, VulcanAction, setStatusFromStatuses, updateFromState} from '../../actions/vulcan_actions';
import {defaultApiHelpers} from '../../contexts/api';
import {VulcanState} from '../../reducers/vulcan_reducer';
import {allFilesToBuffer, filesReturnToMultiFileContent} from '../../selectors/workflow_selectors';
import {MultiFileContent, RunStatus} from '../../api_types';

import {runPromise, useAsyncCallback} from 'etna-js/utils/cancellable_helpers';

const refreshSuggestionText = 'Try Refreshing to fix the previous error.  Alert the Data Library Team if this continues to happen.';

const imageExtRegEx = /\.(png|tiff?|jpe?g|bmp|svg|gif|webp)$/
const largeFileExtRegEx = /\.(gz|zip|Rds|rds|Rdata|RData|rdata|h5ad|h5mu|h5)$/

export function useDataSync(
    state: VulcanState,
    dispatch: Dispatch<VulcanAction>,
    showError: typeof defaultApiHelpers.showError,
    showErrors: typeof defaultApiHelpers.showErrors,
    getFileNames: typeof defaultApiHelpers.getFileNames,
    readFiles: typeof defaultApiHelpers.readFiles,
    postUIValues: typeof defaultApiHelpers.postUIValues,
    getState: typeof defaultApiHelpers.getState,
) {
  const {update_state, configId, update_files, workspace, workspaceId, projectName, status, pushSteps} = state;

  useEffect(() => {
    if (pushSteps.length > 0 && !!workspaceId) {
      // Push files / params to compute server for step
      const pushStep = [...pushSteps][0]
      dispatch(startSyncing())
      showErrors(postUIValues(projectName,workspaceId,status,pushStep), (e) => {
        dispatch(endSyncing(pushStep));
        showError(refreshSuggestionText);
      })
      .then((accountingResponse) => {
        dispatch(useUIAccounting(accountingResponse, pushStep));
      })
    }
  }, [pushSteps, postUIValues])

  useEffect(() => {
    if (update_state) {
      if (configId==null) {
        console.warn('configId unknown, cannot probe state')
        return
      }
      console.log("Checking Current State")
      dispatch(startSyncing())
      showErrors(getState(projectName, configId), (e) => {
        dispatch(endSyncing());
        showError(refreshSuggestionText);
      })
      .then((stateResponse) => {
        dispatch(updateFromState(stateResponse));
      })
    }
  }, [update_state, projectName, configId, refreshSuggestionText, getState])

  useEffect(() => {
    if (update_files && !!workspaceId && !!projectName) {
      console.log("Using data buffering")
      dispatch(startSyncing())
      let update = {output_files: [] as string[], file_contents: {} as {[k: string]: any}};
      showErrors(getFileNames(projectName, workspaceId), (e) => {dispatch(updateFiles(update))})
      .then((fileNamesRaw) => {
        // We want to ignore files that snakemake determined will need to be re-made
        const staleFiles = status.last_file_accounting.planned.concat(status.last_file_accounting.unscheduled);
        const usableFileNames = fileNamesRaw.files.filter(f => !staleFiles.includes(`output/${f}`));
        update['output_files'] = usableFileNames;
        let filesContent: MultiFileContent = {};
          
        const filesReady = allFilesToBuffer(workspace).filter(f => usableFileNames.includes(f));
        // Warn if considering files not reported by snakemake
        const unknownFileNames = filesReady.filter(f => !status.last_file_accounting.completed.includes(`output/${f}`));
        if (unknownFileNames.length>0) {
          console.log(`Note: Utilizing files unknown to snakemake. ${unknownFileNames.join(', ')}`)
        }
        if (filesReady.length == 0) {
          dispatch(updateFiles(update));
        } else {
          // Stub images
          const images = filesReady.filter(f => imageExtRegEx.test(f));
          if (images.length > 0) {
            for (let ind in images) {
              filesContent[images[ind]] = '__IMAGE_CONTENT_STUB__'
            }
          }
          // Stub files meant only for download
          const largeFiles = filesReady.filter(f => largeFileExtRegEx.test(f));
          if (largeFiles.length > 0) {
            for (let ind in largeFiles) {
              filesContent[largeFiles[ind]] = '__LARGE_FILE_CONTENT_STUB__'
            }
          }
          // Grab non-stubbed files
          const files = filesReady.filter(f => !Object.keys(filesContent).includes(f));
          if (files.length > 0) {
            showErrors(readFiles(projectName, workspaceId, files), (e) => {dispatch(updateFiles(update))})
            .then((filesContentRaw) => {
              update['file_contents'] = {
                ...filesContent,
                ...filesReturnToMultiFileContent(filesContentRaw)
              }
              dispatch(updateFiles(update))
            })
          } else {
            dispatch(updateFiles(update))
          }
        }
      })
    }
  }, [update_files, workspaceId, projectName, getFileNames, readFiles]);
}

function delay(ms: number) {
  return new Promise(resolve => setTimeout(resolve, ms));
}

export function useRunSyncing(
  projectName: string,
  workspaceId: number | null,
  runId: number | null,
  showError: typeof defaultApiHelpers.showError,
  showErrors: typeof defaultApiHelpers.showErrors,
  pullRunStatus: typeof defaultApiHelpers.pullRunStatus,
  getIsRunning: typeof defaultApiHelpers.getIsRunning,
  dispatch: Dispatch<VulcanAction>,
): {
  requestRunPolling: () => Promise<unknown>,
  cancelRunning: () => {}
} {
  function suggestRefresh(e: any) {showError(refreshSuggestionText)};
  const [requestRunPolling, cancelRunning] = useAsyncCallback(function* () {
    if (!workspaceId || !runId) {
      showError('Possible UI Bug?: Missing info needed for checking workspace run')
      return
    };
    let isRunning = true;
    while (isRunning) {
      yield delay(3000);
      const runningResponse = yield* runPromise(showErrors(getIsRunning(projectName,workspaceId), suggestRefresh));
      isRunning = runningResponse['running'];
      const response: RunStatus = yield* runPromise(showErrors(pullRunStatus(projectName,workspaceId,runId), suggestRefresh));
      dispatch(setStatusFromStatuses(response, isRunning));
    }
  },
  [projectName, workspaceId, runId, showError, showErrors, pullRunStatus, getIsRunning, dispatch],
  () => {
    // ToDo: cancelRunning!
  });

  return {
    requestRunPolling, cancelRunning
  };
}
