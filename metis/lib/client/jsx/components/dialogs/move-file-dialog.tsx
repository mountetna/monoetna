import React, {useState, useCallback} from 'react';

import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';

import ConfigRow from './config-row';

const MoveFileDialog = ({
  currentBucketName,
  onSubmit
}: {
  currentBucketName: string;
  onSubmit: (bucketName: string, filePath: string) => void;
}) => {
  const [bucketName, setBucketName] = useState(currentBucketName);
  const [newFilePath, setNewFilePath] = useState('');
  const invoke = useActionInvoker();

  const submit = useCallback(() => {
    onSubmit(bucketName, newFilePath);
    invoke({type: 'DISMISS_DIALOG'});
  }, [bucketName, newFilePath]);

  return (
    <div className='move-file-dialog'>
      <div className='title'>Move file</div>
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
          value={newFilePath}
          onChange={(e) => setNewFilePath(e.target.value)}
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

export default MoveFileDialog;
