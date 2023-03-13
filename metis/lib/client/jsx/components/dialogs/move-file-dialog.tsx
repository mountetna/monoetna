import React, {useState, useCallback} from 'react';

import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
import {PickBucket, PickFolder} from 'etna-js/components/metis_exploration';

const MoveFileDialog = ({
  currentBucketName,
  currentPath,
  fileName,
  onSubmit
}: {
  currentBucketName: string;
  currentPath: string;
  fileName: string;
  onSubmit: (bucketName: string, filePath: string) => void;
}) => {
  const [bucketName, setBucketName] = useState(currentBucketName);
  const [newFilePath, setNewFilePath] = useState(currentPath.split("/").slice(0,-1).join("/"));
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
      <div className='title'>Move file: {fileName}</div>
      <PickBucket
        bucket={bucketName}
        label="Bucket"
        setBucket={(e: any) => changeBucket(e)}
      />
      <PickFolder
        bucket={bucketName}
        label="Destination Folder"
        setPath={(e: any) => setNewFilePath(e)}
        basePath={''}
        topLevelPlaceholder={'top-level of bucket'}
        path={newFilePath}
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
