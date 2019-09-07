import * as React from 'react';
import { connect } from 'react-redux';
import { ListEntryColumn, ListEntryTypeColumn, ListEntryUpdatedColumn } from './list-entry';
import { FolderLink } from '../folder-link';
import BucketControl from './bucket-control';
import { selectUserRole } from '../../selectors/user-selector';

const ListBucket = ({bucket, widths, role}) => (
  <div className='list-entry-group'>
    <ListEntryTypeColumn icon='trash' widths={ widths } />
    <ListEntryColumn className='name' widths={widths}>
      <div className='list-entry-file-name' title={bucket.bucket_name}>
        <FolderLink bucket_name={bucket.bucket_name} />
      </div>
    </ListEntryColumn>
    <ListEntryColumn className='description' widths={widths}>
      { bucket.description }
    </ListEntryColumn>
    <ListEntryColumn className='access' widths={widths}>
      <span title={bucket.access}>{ bucket.access }</span>
    </ListEntryColumn>
    <ListEntryColumn className='size' widths={widths}>
      { bucket.count } files
    </ListEntryColumn>
    <ListEntryColumn className='control' widths={widths}>
      {
        role == 'administrator' && <BucketControl bucket={bucket}/>
      }
    </ListEntryColumn>
  </div>
);

export default connect(
  (state) => ({ role: selectUserRole(state) })
)(ListBucket);
