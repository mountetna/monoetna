import * as React from 'react';
import { ListEntryColumn, ListEntryTypeColumn, ListEntryUpdatedColumn } from 'etna-js/components/ListEntry';
import FileControl from './file-control';
import Icon from 'etna-js/components/icon';
import { byteFormat } from '../../utils/format';

const basename = (path) => path.split(/\//).pop();
const foldername = (file) => {
  let [ basename, ...folder_names ] = file.file_name.split(/\//).reverse();
  return folder_names.reverse().join('/');
};

const ListEntryFileNameColumn = ({file, widths}) => (
  <ListEntryColumn className='name' widths={widths}>
    <div className='list-entry-file-name' title={file.file_name}>
      <a href={file.download_url}>{basename(file.file_name)}</a>
    </div>
  </ListEntryColumn>
);

const ListEntryFileStatusColumn = ({file, widths}) => (
  <ListEntryColumn className='status' widths={widths}>
    { file.file_hash && !file.file_hash.match(/^temp-/) && <Icon icon='shield-alt' title={`MD5: ${file.file_hash}`}/> }
    { file.archive_id && <Icon icon='cubes' title='Backed up'/> }
  </ListEntryColumn>
);

const ListEntryFileTypeColumn = ({file: { read_only },widths}) =>
  !read_only
    ? <ListEntryTypeColumn icon='file-alt' widths={widths}/>
    : <ListEntryColumn className='type' widths={widths}>
        <Icon icon='file-alt' overlay='lock'/>
      </ListEntryColumn>;

const ListEntrySizeColumn = ({file,widths}) =>
  <ListEntryColumn className='size' widths={widths}>
    <div className='list-entry-file-size'>
      {byteFormat(file.size, true)}
    </div>
  </ListEntryColumn>;

const ListFile = ({file,current_folder,bucket_name,widths}) => (
  <div className='list-entry-group'>
    <ListEntryFileTypeColumn widths={widths} file={file}/>
    <ListEntryFileNameColumn file={file} widths={widths}/>
    <ListEntryFileStatusColumn file={file} widths={widths}/>
    <ListEntryUpdatedColumn obj={file} widths={widths}/>
    <ListEntrySizeColumn file={file} widths={widths}/>
    <ListEntryColumn className='control' widths={widths}>
      <FileControl file={file} bucket_name={bucket_name} current_folder={current_folder} />
    </ListEntryColumn>
  </div>
);

export default ListFile;
