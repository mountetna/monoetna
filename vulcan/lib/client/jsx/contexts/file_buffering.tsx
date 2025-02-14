import {Dispatch} from 'react';
import {updateFiles, VulcanAction} from '../actions/vulcan_actions';
import {defaultApiHelpers} from './api';
import {VulcanState} from '../reducers/vulcan_reducer';
import {allFilesToBuffer, filesReturnToMultiFileContent} from '../selectors/workflow_selectors';
import {useAsync} from 'etna-js/utils/cancellable_helpers';
import { MultiFileContent } from '../api_types';

const imageExtRegEx = /\.(png|tiff?|jpe?g|bmp|svg|gif|webp)$/

export function useFileBuffering(
    state: VulcanState,
    dispatch: Dispatch<VulcanAction>,
    showErrors: typeof defaultApiHelpers.showErrors,
    getFileNames: typeof defaultApiHelpers.getFileNames,
    readFiles: typeof defaultApiHelpers.readFiles,
) {
  const {update_files, workspace, workspaceId, projectName} = state;

  useAsync(function* () {
    if (update_files && !!workspaceId && !!projectName) {
      console.log("Using data buffering")
      let update = {output_files: [] as string[], file_contents: {} as {[k: string]: any}};
      showErrors(getFileNames(projectName, workspaceId))
      .then((fileNamesRaw) => {
        const fileNames = fileNamesRaw.files;
        update['output_files'] = fileNames;
        let filesContent: MultiFileContent = {};
          
        const filesReady = allFilesToBuffer(workspace).filter(f => fileNames.includes(f));
        if (filesReady.length > 0) {
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
            showErrors(readFiles(projectName, workspaceId, files))
            .then((filesContentRaw) => {
              update['file_contents'] = {
                ...update['file_contents'],
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
  }, [update_files, workspaceId, projectName]);
}
