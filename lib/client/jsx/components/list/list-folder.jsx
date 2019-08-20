import * as React from 'react';
import { ListEntryColumn, ListEntryTypeColumn, ListEntryUpdatedColumn } from './list-entry';
import FolderControl from './folder-control';
import { FolderLink } from '../folder-link';

const ListEntryFolderNameColumn = ({folder, current_bucket, current_folder, widths}) => (
  <ListEntryColumn className='name' widths={widths}>
    <div className='list-entry-file-name' title={folder.folder_name}>
      <FolderLink
        bucket_name={current_bucket}
        folder_path={current_folder}
        folder_name={folder.folder_name}
      />
    </div>
  </ListEntryColumn>
);

const ListFolder = ({ folder, current_folder, current_bucket, widths }) => (
  <div className='list-entry-group'>
    <ListEntryTypeColumn icon='folder' widths={ widths } />
    <ListEntryFolderNameColumn folder={folder}
      current_folder={ current_folder }
      current_bucket={ current_bucket }
      widths={widths} />
    <ListEntryColumn className='status' widths={widths}/>
    <ListEntryUpdatedColumn obj={folder} widths={widths}/>
    <ListEntryColumn className='size' widths={widths}/>
    <ListEntryColumn className='control' widths={widths}>
      <FolderControl folder={folder} current_folder={current_folder} />
    </ListEntryColumn>
  </div>
);

export default ListFolder;
