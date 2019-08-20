import * as React from 'react';
import { ListEntryColumn, ListEntryTypeColumn, ListEntryUpdatedColumn } from './list-entry';
import { FolderLink } from '../folder-link';

const ListBucket = ({bucket, widths}) => (
  <div className='list-entry-group'>
    <ListEntryTypeColumn icon='trash' widths={ widths } />
    <ListEntryColumn className='name' widths={widths}>
      <div className='list-entry-file-name' title={bucket.bucket_name}>
        <FolderLink bucket_name={bucket.bucket_name} />
      </div>
    </ListEntryColumn>
    <ListEntryColumn className='description' widths={widths}>
      { bucket.count } files
    </ListEntryColumn>
  </div>
);

export default ListBucket;
