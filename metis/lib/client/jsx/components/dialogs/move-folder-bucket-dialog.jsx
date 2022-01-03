import React, {useState, useCallback} from 'react';

import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';

import ConfigRow from './config-row';

const MoveFolderBucketDialog = ({onSubmit}) => {
  const [bucketName, setBucketName] = useState(null);
  const [newFolderPath, setNewFolderPath] = useState(null);
  const invoke = useActionInvoker();

  const submit = useCallback(() => {
    onSubmit(bucketName, newFolderPath);
    invoke({type: 'DISMISS_DIALOG'});
  }, [bucketName, newFolderPath]);

  return (
    <div className='move-folder-bucket-dialog'>
      <div className='title'>Move folder to another bucket</div>
      <ConfigRow label='Bucket name'>
        <input
          type='text'
          placeholder='E.g. bucket_name'
          value={bucketName}
          onChange={(e) => setBucketName(e.target.value)}
        />
      </ConfigRow>
      <ConfigRow label='Folder path'>
        <input
          type='text'
          placeholder='E.g. root/child/folder_name'
          value={newFolderPath}
          onChange={(e) => setNewFolderPath(e.target.value)}
        />
      </ConfigRow>
      <div className='submit'>
        <span
          className='button'
          disabled={!(bucketName || newFolderPath)}
          onClick={submit}
        >
          Move
        </span>
      </div>
    </div>
  );
};

export default MoveFolderBucketDialog;
