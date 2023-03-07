import React, {useState, useCallback} from 'react';

import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
import {PickBucket, PickFileOrFolder} from 'etna-js/components/metis_exploration';

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

  const changeBucket = useCallback((e) => {
    setParentFolderPath('')
    setBucketName(e)
  }, [])

  return (
    <div className='copy-folder-dialog'>
      <div className='title'>Copy folder</div>
      <PickBucket
        bucket={bucketName}
        label="Bucket"
        setBucket={(e: any) => changeBucket(e)}
      />
      <PickFileOrFolder
        bucket={bucketName}
        label="Destination Folder"
        setPath={(e: any) => setParentFolderPath(e)}
        basePath={''}
        topLevelPlaceholer={'top-level of bucket'}
        path={parentFolderPath}
        allowFiles={false}
      />
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
