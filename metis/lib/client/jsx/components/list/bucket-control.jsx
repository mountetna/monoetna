import React, {useCallback} from 'react';
import MenuControl from '../menu-control';
import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';

const BucketControl = ({bucket}) => {
  const invoke = useActionInvoker();

  const updateBucket = useCallback(
    (bucket) => {
      invoke({
        type: 'UPDATE_BUCKET',
        bucket
      });
    },
    [invoke]
  );

  const configureBucket = useCallback(() => {
    let {bucket_name, description, access} = bucket;

    let dialog = {
      type: 'configure-bucket',
      updateBucket,
      bucket_name,
      description,
      access
    };
    invoke({
      type: 'SHOW_DIALOG',
      dialog
    });
  }, [invoke, bucket, updateBucket]);

  const destroyBucket = useCallback(() => {
    invoke({
      type: 'DESTROY_BUCKET',
      bucket
    });
  }, [invoke, bucket]);

  let items = [
    {
      label: 'Configure bucket',
      callback: configureBucket,
      show: true,
      role: 'administrator'
    },
    bucket.count == 0 && {
      label: 'Remove bucket',
      callback: destroyBucket,
      role: 'administrator'
    }
  ].filter((_) => _);

  return <MenuControl items={items} />;
};

export default BucketControl;
