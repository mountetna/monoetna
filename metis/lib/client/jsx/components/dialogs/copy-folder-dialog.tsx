import React, {useState, useCallback} from 'react';

import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';

import ConfigRow from './config-row';

const CopyFolderDialog = ({
  currentBucketName,
  onSubmit
}: {
  currentBucketName: string;
  onSubmit: (bucketName: string, folderPath: string) => void;
}) => {
  const [bucketName, setBucketName] = useState(currentBucketName);
  const [parentFolderPath, setParentFolderPath] = useState('');
  const invoke = useActionInvoker();

  const submit = useCallback(() => {
    onSubmit(bucketName, parentFolderPath);
    invoke({type: 'DISMISS_DIALOG'});
  }, [bucketName, parentFolderPath, invoke, onSubmit]);

  return (
    <div className='copy-folder-dialog'>
      <div className='title'>Copy folder</div>
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
          value={parentFolderPath}
          onChange={(e) => setParentFolderPath(e.target.value)}
        />
      </ConfigRow>
      <div className='submit'>
        <span
          className={`button ${bucketName ? '' : 'disabled'}`}
          onClick={submit}
        >
          Copy
        </span>
      </div>
    </div>
  );
};

export default CopyFolderDialog;
