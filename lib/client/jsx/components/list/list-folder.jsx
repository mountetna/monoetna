import * as React from 'react';
import { ListEntryColumn, ListEntryTypeColumn, ListEntryUpdatedColumn } from './list-entry';
import FolderControl from './folder-control';
import Icon from '../icon';
import { FolderLink } from '../folder-link';

const ListEntryFolderTypeColumn = ({folder: { read_only },widths}) =>
  !read_only
    ? <ListEntryTypeColumn icon='folder' widths={widths}/>
    : <ListEntryColumn className='type' widths={widths}>
        <Icon icon='folder' overlay='lock'/>
      </ListEntryColumn>;

const ListEntryFolderNameColumn = ({folder, bucket_name, current_folder, widths}) => (
  <ListEntryColumn className='name' widths={widths}>
    <div className='list-entry-file-name' title={folder.folder_name}>
      <FolderLink
        bucket_name={bucket_name}
        folder_path={current_folder}
        folder_name={folder.folder_name}
      />
    </div>
  </ListEntryColumn>
);

const ListFolder = ({ folder, current_folder, bucket_name, widths }) => (
  <div className='list-entry-group'>
    <ListEntryFolderTypeColumn widths={widths} folder={folder}/>
    <ListEntryFolderNameColumn folder={folder}
      current_folder={ current_folder }
      bucket_name={ bucket_name }
      widths={widths} />
    <ListEntryColumn className='status' widths={widths}/>
    <ListEntryUpdatedColumn obj={folder} widths={widths}/>
    <ListEntryColumn className='size' widths={widths}/>
    <ListEntryColumn className='control' widths={widths}>
      <FolderControl bucket_name={ bucket_name } folder={folder} current_folder={current_folder} />
    </ListEntryColumn>
  </div>
);

export default ListFolder;
