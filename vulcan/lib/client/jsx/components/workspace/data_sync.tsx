import {Dispatch, useEffect} from 'react';
import {removeSync, updateFiles, useUIAccounting, VulcanAction} from '../../actions/vulcan_actions';
import {defaultApiHelpers} from '../../contexts/api';
import {VulcanState} from '../../reducers/vulcan_reducer';
import {allFilesToBuffer, filesReturnToMultiFileContent} from '../../selectors/workflow_selectors';
import { MultiFileContent } from '../../api_types';

const imageExtRegEx = /\.(png|tiff?|jpe?g|bmp|svg|gif|webp)$/

export function useDataSync(
    state: VulcanState,
    dispatch: Dispatch<VulcanAction>,
    showErrors: typeof defaultApiHelpers.showErrors,
    getFileNames: typeof defaultApiHelpers.getFileNames,
    readFiles: typeof defaultApiHelpers.readFiles,
    postUIValues: typeof defaultApiHelpers.postUIValues,
) {
  const {update_files, workspace, workspaceId, projectName, status, pushSteps} = state;

  useEffect(() => {
    if (pushSteps.length > 0 && !!workspaceId) {
      // Push files / params to compute server for step
      const pushStep = [...pushSteps][0]
      showErrors(postUIValues(projectName,workspaceId,status,pushStep))
      .then((accountingResponse) => {
        dispatch(useUIAccounting(accountingResponse, pushStep));
        dispatch(removeSync(pushStep));
      })
    }
  }, [pushSteps, postUIValues])

  useEffect(() => {
    if (update_files && !!workspaceId && !!projectName) {
      console.log("Using data buffering")
      let update = {output_files: [] as string[], file_contents: {} as {[k: string]: any}};
      showErrors(getFileNames(projectName, workspaceId), (e) => {dispatch(updateFiles(update))})
      .then((fileNamesRaw) => {
        const fileNames = fileNamesRaw.files;
        update['output_files'] = fileNames;
        let filesContent: MultiFileContent = {};
          
        const filesReady = allFilesToBuffer(workspace).filter(f => fileNames.includes(f));
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
