import React, {useState, useCallback} from 'react';

import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
import {PickBucket, PickFolder} from 'etna-js/components/metis_exploration';

const CopyFolderDialog = ({
  currentBucketName,
  currentPath,
  folderName,
  onSubmit
}: {
  currentBucketName: string;
  currentPath: string;
  folderName: string;
  onSubmit: (bucketName: string, folderPath: string) => void;
}) => {
  const [bucketName, setBucketName] = useState(currentBucketName);
  const [parentFolderPath, setParentFolderPath] = useState(currentPath.split("/").slice(0,-1).join("/"));
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
      <div className='title'>Copy folder: {folderName}</div>
      <PickBucket
        bucket={bucketName}
        label="Bucket"
        setBucket={(e: any) => changeBucket(e)}
      />
      <PickFolder
        bucket={bucketName}
        label="Destination Folder"
        setPath={(e: any) => setParentFolderPath(e)}
        basePath={''}
        topLevelPlaceholder={'top-level of bucket'}
        path={parentFolderPath}
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
