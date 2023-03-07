import React, {useState, useCallback} from 'react';

import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
import {PickBucket, PickFileOrFolder} from 'etna-js/components/metis_exploration';

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
  }, [bucketName, newFilePath, invoke, onSubmit]);

  const changeBucket = useCallback((e) => {
    setNewFilePath('')
    setBucketName(e)
  }, [])

  return (
    <div className='move-file-dialog'>
      <div className='title'>Move file</div>
      <PickBucket
        bucket={bucketName}
        label="Bucket"
        setBucket={(e: any) => changeBucket(e)}
      />
      <PickFileOrFolder
        bucket={bucketName}
        label="Destination Folder"
        setPath={(e: any) => setNewFilePath(e)}
        basePath={''}
        topLevelPlaceholer={'top-level of bucket'}
        path={newFilePath}
        allowFiles={false}
      />
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
