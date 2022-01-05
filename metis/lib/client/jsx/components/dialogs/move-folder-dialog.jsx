import React, {useState, useCallback} from 'react';

import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';

import ConfigRow from './config-row';

const MoveFolderDialog = ({currentBucketName, onSubmit}) => {
  const [bucketName, setBucketName] = useState(currentBucketName);
  const [newFolderPath, setNewFolderPath] = useState('');
  const invoke = useActionInvoker();

  const submit = useCallback(() => {
    onSubmit(bucketName, newFolderPath);
    invoke({type: 'DISMISS_DIALOG'});
  }, [bucketName, newFolderPath]);

  return (
    <div className='move-folder-dialog'>
      <div className='title'>Move folder</div>
      <ConfigRow label='Bucket name'>
        <input
          type='text'
          placeholder='E.g. bucket_name'
          value={bucketName}
          onChange={(e) => setBucketName(e.target.value)}
        />
      </ConfigRow>
      <ConfigRow label='Parent folder (blank for root)'>
        <input
          type='text'
          placeholder='E.g. root/child'
          value={newFolderPath}
          onChange={(e) => setNewFolderPath(e.target.value)}
        />
      </ConfigRow>
      <div className='submit'>
        <span className='button' disabled={!bucketName} onClick={submit}>
          Move
        </span>
      </div>
    </div>
  );
};

export default MoveFolderDialog;
