import React, {useCallback} from 'react';
import {filePath, includesFolders} from 'etna-js/utils/file';
import MenuControl from '../menu-control';

import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
import {message} from '../../actions/message_actions';
import {Folder, UiControlItem} from '../../types/metis-types';

const FolderControl = ({
  folder,
  current_folder,
  bucket_name
}: {
  folder: Folder;
  current_folder: string;
  bucket_name: string;
}) => {
  const invoke = useActionInvoker();

  const unprotectFolder = useCallback(() => {
    invoke({type: 'UNPROTECT_FOLDER', folder, bucket_name});
  }, [folder, bucket_name]);

  const protectFolder = useCallback(() => {
    invoke({type: 'PROTECT_FOLDER', folder, bucket_name});
  }, [folder, bucket_name]);

  const renameFolder = useCallback(() => {
    let new_folder_name = prompt(
      'What is the new name of this folder?',
      folder.folder_name
    );
    if (includesFolders(new_folder_name)) {
      invoke(
        message(
          'warning',
          'Folder renaming failed',
          'Invalid name -- cannot change the folder path'
        )
      );
    } else if (new_folder_name)
      invoke({
        type: 'RENAME_FOLDER',
        folder,
        new_folder_path: filePath(current_folder, new_folder_name),
        bucket_name
      });
  }, [folder, current_folder, bucket_name]);

  const removeFolder = useCallback(() => {
    invoke({type: 'REMOVE_FOLDER', folder, bucket_name});
  }, [folder, bucket_name]);

  const touchFolder = useCallback(() => {
    invoke({type: 'TOUCH_FOLDER', folder, bucket_name});
  }, [folder, bucket_name]);

  const moveFolder = useCallback(
    (newBucketName, newFolderPath) => {
      invoke({
        type: 'RENAME_FOLDER',
        folder,
        new_folder_path: filePath(newFolderPath, folder.folder_name),
        new_bucket_name: newBucketName,
        bucket_name,
        current_folder: current_folder || ''
      });
    },
    [folder, bucket_name, current_folder]
  );

  const moveFolderDialog = useCallback(() => {
    let dialog = {
      type: 'move-folder',
      onSubmit: moveFolder,
      currentBucketName: bucket_name
    };
    invoke({
      type: 'SHOW_DIALOG',
      dialog
    });
  }, [moveFolder, bucket_name]);

  let items: UiControlItem[] = folder.read_only
    ? [
        {
          label: 'Unprotect folder',
          callback: unprotectFolder,
          role: 'administrator'
        }
      ]
    : [
        {
          label: 'Rename folder',
          callback: renameFolder,
          role: 'editor'
        },
        {
          label: 'Protect folder',
          callback: protectFolder,
          role: 'administrator'
        },
        {
          label: 'Remove folder',
          callback: removeFolder,
          role: 'editor'
        },
        {
          label: 'Touch folder',
          callback: touchFolder,
          role: 'editor'
        },
        {
          label: 'Move folder',
          callback: moveFolderDialog,
          role: 'editor'
        }
      ];
  return <MenuControl items={items} />;
};

export default FolderControl;
