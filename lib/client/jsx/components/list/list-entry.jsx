import * as React from 'react';
import { userFormat, byteFormat, dateFormat } from '../../utils/format';

import UploadMeter from './upload-meter';
import UploadControl from './upload-control';
import FileControl from './file-control';
import { FolderLink } from '../folder-link';
import Icon from '../icon';
import { filePath } from '../../utils/file';

const ListEntryColumn = ({className,widths,children}) =>
  <div className={`list-entry-column-group ${className}`}
    style={ { flexBasis: widths[className] } } >
    { children }
  </div>;

const ListEntryTypeColumn = ({icon, widths}) =>
  <ListEntryColumn className='type' widths={widths}>
    <Icon icon={icon}/>
  </ListEntryColumn>;

const ListEntryControlColumn = ({widths, ...props}) =>
  <ListEntryColumn className='control' widths={widths}>
    <FileControl { ...props } />
  </ListEntryColumn>;

const ListEntrySizeColumn = ({file,widths}) =>
  <ListEntryColumn className='size' widths={widths}>
    <div className='list-entry-file-size'>
      {byteFormat(file.size, true)}
    </div>
  </ListEntryColumn>;

const ListEntryUpdatedColumn = ({file, widths}) =>
  <ListEntryColumn className='updated' widths={widths}>
    <div className='list-entry-updated-name'>
      {dateFormat(file.updated_at)} by {userFormat(file.author)}
    </div>
    <div className='list-entry-role'>
      {file.role}
    </div>
  </ListEntryColumn>;

const basename = (path) => path.split(/\//).pop();
const foldername = (file) => {
  let [ basename, ...folder_names ] = file.file_name.split(/\//).reverse();
  return folder_names.reverse().join('/');
}

const browse_path = (folder_name, current_folder) => (
  `/${CONFIG.project_name}/browse/${filePath(current_folder, folder_name)}`
);

const ListEntryFoldernameColumn = ({folder, current_folder, widths}) => (
  <ListEntryColumn className='name' widths={widths}>
    <div className='list-entry-file-name' title={folder.folder_name}>
      <FolderLink
        folder_name={folder.folder_name}
        folder_path={current_folder}/>
    </div>
  </ListEntryColumn>
)

const ListEntryFileNameColumn = ({file, widths}) => (
  <ListEntryColumn className='name' widths={widths}>
    <div className='list-entry-file-name' title={file.file_name}>
      <a href={file.download_url}>{basename(file.file_name)}</a>
    </div>
  </ListEntryColumn>
)

const ListEntryFileTypeColumn = ({file,widths}) => {
  let { read_only } = file;

  if (!read_only)
    return <ListEntryTypeColumn icon='file-alt' widths={widths}/>;

  return <ListEntryColumn className='type' widths={widths}>
    <Icon icon='file-alt' overlay='lock'/>
  </ListEntryColumn>;

};

const ListEntryFileStatusColumn = ({file, widths}) => (
  <ListEntryColumn className='status' widths={widths}>
    { file.file_hash && <Icon icon='shield-alt' title={`MD5: ${file.file_hash}`}/> }
    { file.archive_id && <Icon icon='cubes' title='Backed up'/> }
  </ListEntryColumn>
)


export const ListFile = ({file,current_folder,widths}) => (
  <div className='list-entry-group'>
    <ListEntryFileTypeColumn widths={widths} file={file}/>
    <ListEntryFileNameColumn file={file} widths={widths}/>
    <ListEntryFileStatusColumn file={file} widths={widths}/>
    <ListEntryUpdatedColumn file={file} widths={widths}/>
    <ListEntrySizeColumn file={file} widths={widths}/>
    <ListEntryControlColumn file={file} current_folder={current_folder} widths={widths}/>
  </div>
);

export const ListFolder = ({ folder, current_folder, widths }) => (
  <div className='list-entry-group'>
    <ListEntryTypeColumn icon='folder' widths={ widths } />
    <ListEntryFoldernameColumn folder={folder}
      current_folder={ current_folder } widths={widths} />
    <ListEntryColumn className='status' widths={widths}/>
    <ListEntryUpdatedColumn file={folder} widths={widths}/>
    <ListEntryColumn className='size' widths={widths}/>
    <ListEntryColumn className='control' widths={widths}/>
  </div>
);

export const ListUpload = ({ upload, widths }) => (
  <div className='list-entry-group upload'>
    <ListEntryTypeColumn icon='upload' widths={ widths } />
    <ListEntryColumn className='name' widths={ widths }>
      <span className='list-entry-file-name'>
        {upload.file_name}
      </span>
    </ListEntryColumn>
    <ListEntryColumn className='status' widths={widths}/>
    <ListEntryColumn className='updated' widths={widths}>
      <UploadMeter upload={ upload }/>
    </ListEntryColumn>
    <ListEntryColumn className='size' widths={widths}/>
    <ListEntryColumn className='control' widths={widths}>
      <UploadControl upload={ upload } />
    </ListEntryColumn>
  </div>
);
