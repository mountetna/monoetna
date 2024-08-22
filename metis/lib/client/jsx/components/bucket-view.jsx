import React, {useEffect, useState} from 'react';

import {useReduxState} from 'etna-js/hooks/useReduxState';
import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
import {useAsync} from 'etna-js/utils/cancellable_helpers'

import ListBucket from './list/list-bucket';
import ListHead from './list/list-head';
import FolderBreadcrumb from './folder-breadcrumb';
import ControlBar from './control-bar';

import {selectBuckets} from 'etna-js/selectors/directory-selector';

import AutorenewIcon from '@material-ui/icons/Autorenew';
import { makeStyles } from '@material-ui/core/styles';

const useStyles = makeStyles((theme) => ({
  loadingIcon: {
    'animation': '$spin 4s linear infinite'
  },
  '@keyframes spin': {
      '100%': {
          '-webkit-transform': 'rotate(360deg)',
          'transform': 'rotate(360deg)',
      }
  },
}))

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

  const [loading, setLoading] = useState(true)
  const classes = useStyles();

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

  useAsync( async () => {
    await invoke({type: 'RETRIEVE_BUCKETS'});
    setLoading(false);
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
        {loading ? <>
          <AutorenewIcon className={classes.loadingIcon}/>
          Loading
        </> : null}
      </div>
    </div>
  );
};

export default BucketView;
