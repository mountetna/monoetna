import React, {useCallback} from 'react';
import {filePath} from 'etna-js/utils/file';
import MenuControl from '../menu-control';

import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';

const FolderControl = ({folder, current_folder, bucket_name}) => {
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
    if (new_folder_name)
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

  let items = folder.read_only
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
        }
      ];
  return <MenuControl items={items} />;
};

export default FolderControl;
