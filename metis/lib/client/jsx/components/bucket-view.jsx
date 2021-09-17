import React, {useEffect} from 'react';

import {useReduxState} from 'etna-js/hooks/useReduxState';
import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';

import ListBucket from './list/list-bucket';
import ListHead from './list/list-head';
import FolderBreadcrumb from './folder-breadcrumb';
import ControlBar from './control-bar';

import {selectBuckets} from 'etna-js/selectors/directory-selector';

const COLUMNS = [
  {name: 'type', width: '60px'},
  {name: 'name', width: '30%'},
  {name: 'description', width: '50%'},
  {
    name: 'access',
    width: '10%',
    title: 'Permission level or access list required to use this bucket'
  },
  {name: 'size', width: '10%'},
  {name: 'control', width: '100px', hide: true}
];

const COLUMN_WIDTHS = COLUMNS.reduce((widths, column) => {
  widths[column.name] = column.width;
  return widths;
}, {});

const BucketView = () => {
  const invoke = useActionInvoker();
  const buckets = useReduxState((state) => selectBuckets(state));

  function createBucket(bucket) {
    invoke({type: 'CREATE_BUCKET', bucket});
  }

  function onShowDialog() {
    let dialog = {
      type: 'configure-bucket',
      createBucket
    };
    invoke({type: 'SHOW_DIALOG', dialog});
  }

  let buttons = [
    {
      onClick: onShowDialog,
      title: 'Create bucket',
      icon: 'trash',
      overlay: 'plus',
      role: 'administrator'
    }
  ];

  useEffect(() => {
    invoke({type: 'RETRIEVE_BUCKETS'});
  }, []);

  return (
    <div className='bucket-view-group'>
      <div className='control-group'>
        <FolderBreadcrumb />
        <ControlBar buttons={buttons} />
      </div>
      <ListHead columns={COLUMNS} />
      <div id='list-body-group'>
        {Object.values(buckets)
          .sort((b) => b.created_at)
          .map((bucket) => (
            <ListBucket
              key={bucket.bucket_name}
              widths={COLUMN_WIDTHS}
              bucket={bucket}
            />
          ))}
      </div>
    </div>
  );
};

export default BucketView;
