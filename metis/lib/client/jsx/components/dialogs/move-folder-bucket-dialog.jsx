import React, {useCallback} from 'react';

const MoveFolderBucketDialog = ({dismissDialog}) => {
  return (
    <div className='move-folder-bucket-dialog'>
      <div className='title'>Move folder to another bucket</div>
      <div className=''>
        <label>
          Bucket name
          <input type='text' name='bucket-name' />
        </label>
        <label>
          Folder path
          <input type='text' name='folder-path' />
        </label>
      </div>
    </div>
  );
};

export default MoveFolderBucketDialog;
