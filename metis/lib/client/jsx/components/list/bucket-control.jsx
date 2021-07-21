import React, {useCallback} from 'react';
import MenuControl from '../menu-control';
import {selectUser} from 'etna-js/selectors/user-selector';
import {useReduxState} from 'etna-js/hooks/useReduxState';
import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
import {isAdmin} from 'etna-js/utils/janus';
import {useFeatureFlag} from 'etna-js/hooks/useFeatureFlag';

const BucketControl = ({bucket}) => {
  const user = useReduxState((state) => selectUser(state));
  const canIngest = useFeatureFlag('ingest');
  const invoke = useActionInvoker();

  function updateBucket(bucket) {
    invoke({
      type: 'UPDATE_BUCKET',
      bucket
    });
  }

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
  }, [bucket.bucket_name, bucket.description, bucket.access, invoke]);

  const destroyBucket = useCallback(() => {
    invoke({
      type: 'DESTROY_BUCKET',
      bucket
    });
  }, [invoke, bucket]);

  const showIngestDialog = useCallback(() => {
    let dialog = {
      type: 'ingest-to-bucket'
    };
    invoke({
      type: 'SHOW_DIALOG',
      dialog
    });
  }, [invoke]);

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

  if (isAdmin(user, CONFIG.project_name) && canIngest) {
    items.push({
      label: 'Ingest files',
      callback: showIngestDialog
    });
  }

  return <MenuControl items={items} />;
};

export default BucketControl;
