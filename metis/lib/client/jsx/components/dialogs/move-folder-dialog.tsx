import React, {useState, useCallback} from 'react';

import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
import {PickBucket} from 'etna-js/components/metis_exploration';

import ConfigRow from './config-row';

const MoveFolderDialog = ({
  currentBucketName,
  onSubmit
}: {
  currentBucketName: string;
  onSubmit: (bucketName: string, folderPath: string) => void;
}) => {
  const [bucketName, setBucketName] = useState(currentBucketName);
  const [newFolderPath, setNewFolderPath] = useState('');
  const invoke = useActionInvoker();

  const submit = useCallback(() => {
    onSubmit(bucketName, newFolderPath);
    invoke({type: 'DISMISS_DIALOG'});
  }, [bucketName, newFolderPath, invoke, onSubmit]);

  return (
    <div className='move-folder-dialog'>
      <div className='title'>Move folder</div>
      <ConfigRow label='Bucket:'>
        <PickBucket
          setBucket={(e) => setBucketName(e)}
          initialValue={currentBucketName}
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
        <span
          className={`button ${bucketName ? '' : 'disabled'}`}
          onClick={submit}
        >
          Move
        </span>
      </div>
    </div>
  );
};

export default MoveFolderDialog;
