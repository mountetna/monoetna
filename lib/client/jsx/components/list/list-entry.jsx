import * as React from 'react';
import { userFormat, byteFormat, dateFormat } from '../../utils/format';

import UploadMeter from './upload-meter';
import UploadControl from './upload-control';
import FileControl from './file-control';

const ListEntryColumn = ({className,widths,children}) =>
  <div className={`list-entry-column-group ${className}`}
    style={ { flexBasis: widths[className] } } >
    { children }
  </div>;

const ListEntryTypeColumn = ({icon, widths}) =>
  <ListEntryColumn className='type' widths={widths}>
    <span className={ icon
    }/>
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

const basename = (file) => file.file_name.split(/\//).pop();
const foldername = (file) => {
  let [ basename, ...folder_names ] = file.file_name.split(/\//).reverse();
  return folder_names.reverse().join('/');
}
const FolderLink = ({className, file}) =>
    <div className={className} title={file.file_name}>
      <a href={
        `/${CONFIG.project_name}/browse/${file.file_name}`
      }>{basename(file)}</a>
    </div>;

const FileLink = ({className, file}) =>
    <div className={className} title={file.file_name}>
      <a href={file.download_url}>{basename(file)}</a>
    </div>;

const ListEntryFilenameColumn = ({file, widths}) => {
  let LinkType = file.is_folder ? FolderLink : FileLink;
  return <ListEntryColumn className='name' widths={widths}>
    <LinkType className='list-entry-file-name' file={file}/>
  </ListEntryColumn>;
}

export class ListEntry extends React.Component{
  render(){
    let { file, widths } = this.props;

    return (
      <div className='list-entry-group'>
        <ListEntryTypeColumn icon='far fa-file-alt' widths={widths}/>
        <ListEntryFilenameColumn file={file} widths={widths}/>
        <ListEntryUpdatedColumn file={file} widths={widths}/>
        <ListEntrySizeColumn file={file} widths={widths}/>
        <ListEntryControlColumn file={file} widths={widths}/>
      </div>
    );
  }
}

export class ListFolder extends React.Component{
  render(){
    let { file, widths }  = this.props;

    return (
      <div className='list-entry-group'>
        <ListEntryTypeColumn icon='fas fa-folder' widths={ widths } />
        <ListEntryFilenameColumn file={file} widths={widths} />
        <ListEntryUpdatedColumn file={file} widths={widths}/>
        <ListEntryColumn className='size' widths={widths}/>
        <ListEntryColumn className='control' widths={widths}/>
      </div>
    );
  }
}

export class ListUpload extends React.Component{
  render() {
    let { upload, widths } = this.props;

    return (
      <div className='list-entry-group upload'>
        <ListEntryTypeColumn icon='upload fas fa-upload' widths={ widths } />
        <ListEntryColumn className='name' widths={ widths }>
          <span className='list-entry-file-name'>
            {upload.file_name}
          </span>
        </ListEntryColumn>
        <ListEntryColumn className='updated' widths={widths}>
          <UploadMeter upload={ upload }/>
        </ListEntryColumn>
        <ListEntryColumn className='size' widths={widths}/>
        <ListEntryColumn className='control' widths={widths}>
          <UploadControl upload={ upload } />
        </ListEntryColumn>
      </div>
    );
  }
}
