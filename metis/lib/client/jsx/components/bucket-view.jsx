import React, {useContext, useEffect} from 'react';

import {selectUser} from 'etna-js/selectors/user-selector';
import {useReduxState} from 'etna-js/hooks/useReduxState';
import {isAdmin} from 'etna-js/utils/janus';
import {useFeatureFlag} from 'etna-js/hooks/useFeatureFlag';
import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
import {showMessages} from 'etna-js/actions/message_actions';

import ListBucket from './list/list-bucket';
import ListHead from './list/list-head';
import FolderBreadcrumb from './folder-breadcrumb';
import ControlBar from './control-bar';

import {selectBuckets} from 'etna-js/selectors/directory-selector';
import {listHosts} from '../api/polyphemus/ingest_api';
import {IngestContext} from '../contexts/ingest_context';

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
  const user = useReduxState((state) => selectUser(state));
  const canIngest = useFeatureFlag('ingest');
  const {state, setHosts} = useContext(IngestContext);

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

  useEffect(() => {
    const isAdminUser = isAdmin(user, CONFIG.project_name);
    if (
      CONFIG.project_name &&
      isAdminUser &&
      canIngest &&
      0 === Object.keys(state.hosts).length
    ) {
      listHosts()
        .then(({hosts}) => {
          setHosts(hosts);
        })
        .catch((e) => {
          console.error(e);
          invoke(showMessages(e.errors || [e.error] || [e]));
        });
    }
  }, [CONFIG.project_name, user, canIngest, state.hosts, setHosts]);

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
