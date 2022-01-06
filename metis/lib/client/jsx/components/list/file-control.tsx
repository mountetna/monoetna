import React, {useCallback} from 'react';
import {copyText} from 'etna-js/utils/copy';
import {filePath, includesFolders} from 'etna-js/utils/file';
import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
import {message} from '../../actions/message_actions';
import MenuControl from '../menu-control';
import {File, UiControlItem} from '../../types/metis-types';

const FileControl = ({
  file,
  bucket_name,
  current_folder
}: {
  file: File;
  bucket_name: string;
  current_folder: string;
}) => {
  const invoke = useActionInvoker();

  const unprotectFile = useCallback(() => {
    invoke({type: 'UNPROTECT_FILE', file, bucket_name});
  }, [file, bucket_name]);

  const protectFile = useCallback(() => {
    invoke({type: 'PROTECT_FILE', file, bucket_name});
  }, [file, bucket_name]);

  const renameFile = useCallback(() => {
    let new_file_name = prompt(
      'What is the new name of this file?',
      file.file_name
    );
    if (includesFolders(new_file_name)) {
      invoke(
        message(
          'warning',
          'File renaming failed',
          'Invalid name -- cannot change the path'
        )
      );
    } else if (new_file_name)
      invoke({
        type: 'RENAME_FILE',
        file,
        new_file_path: filePath(current_folder, new_file_name),
        bucket_name
      });
  }, [file, current_folder, bucket_name]);

  const removeFile = useCallback(() => {
    invoke({type: 'REMOVE_FILE', file, bucket_name});
  }, [file, bucket_name]);

  const copyLink = useCallback(() => {
    copyText(file.download_url);
  }, [file.download_url, bucket_name]);

  const copyMetisPath = useCallback(() => {
    copyText(
      `metis://${file.project_name}/${file.bucket_name}/${file.file_path}`
    );
  }, [file, bucket_name]);

  const touchFile = useCallback(() => {
    invoke({type: 'TOUCH_FILE', file, bucket_name});
  }, [file, bucket_name]);

  const downloadFile = useCallback(() => {
    let download = document.createElement('a');
    download.setAttribute('href', file.download_url);
    download.setAttribute('download', file.file_name);

    download.style.display = 'none';
    document.body.appendChild(download);

    download.click();

    document.body.removeChild(download);
  }, [file]);

  const moveFile = useCallback(
    (newBucketName, newFilePath) => {
      invoke({
        type: 'RENAME_FILE',
        file,
        new_file_path: filePath(newFilePath, file.file_name),
        new_bucket_name: newBucketName,
        bucket_name,
        current_folder: current_folder || ''
      });
    },
    [file, bucket_name, current_folder]
  );

  const moveFileDialog = useCallback(() => {
    let dialog = {
      type: 'move-file',
      onSubmit: moveFile,
      currentBucketName: bucket_name
    };
    invoke({
      type: 'SHOW_DIALOG',
      dialog
    });
  }, [moveFile, bucket_name]);

  let items: UiControlItem[] = [
    {label: 'Download file', callback: downloadFile, role: 'viewer'},
    {label: 'Copy download link', callback: copyLink, role: 'viewer'},
    {label: 'Copy metis path', callback: copyMetisPath, role: 'viewer'}
  ].concat(
    file.read_only
      ? [
          {
            label: 'Unprotect file',
            callback: unprotectFile,
            role: 'administrator'
          }
        ]
      : [
          {
            label: 'Rename file',
            callback: renameFile,
            role: 'editor'
          },
          {
            label: 'Protect file',
            callback: protectFile,
            role: 'administrator'
          },
          {
            label: 'Remove file',
            callback: removeFile,
            role: 'editor'
          },
          {
            label: 'Touch file',
            callback: touchFile,
            role: 'editor'
          },
          {
            label: 'Move file',
            callback: moveFileDialog,
            role: 'editor'
          }
        ]
  );
  return <MenuControl items={items} />;
};

export default FileControl;
