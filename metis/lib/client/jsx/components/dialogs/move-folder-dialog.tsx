import React, {useState, useCallback} from 'react';

import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
import {PickBucket, PickFileOrFolder} from 'etna-js/components/metis_exploration';

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

  const changeBucket = useCallback((e) => {
    setNewFolderPath('')
    setBucketName(e)
  }, [])

  return (
    <div className='move-folder-dialog'>
      <div className='title'>Move folder</div>
      <PickBucket
        bucket={bucketName}
        label="Bucket"
        setBucket={(e: any) => changeBucket(e)}
      />
      <PickFileOrFolder
        bucket={bucketName}
        label="Destination Folder"
        setPath={(e: any) => setNewFolderPath(e)}
        basePath={''}
        topLevelPlaceholer={'top-level of bucket'}
        path={newFolderPath}
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

export default MoveFolderDialog;
