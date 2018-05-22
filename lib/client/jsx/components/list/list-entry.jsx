import * as React from 'react';
import { userFormat, byteFormat, dateFormat } from '../../utils/format';

import FileControl from './file-control';

const ListEntryColumn = ({className,widths,children}) =>
  <div className={`list-entry-column-group ${className}`}
    style={ { flexBasis: widths[className] } } >
    { children }
  </div>;

const ListEntryTypeColumn = ({is_folder, widths}) =>
  <ListEntryColumn className='type' widths={widths}>
    <span className={ is_folder ? 'fas fa-folder' : 'far fa-file-alt' }/>
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
  constructor(){
    super();
  }

  removeFile(){
    this.props.callbacks.removeFile(this.props.file);
  }

  render(){
    let { file, widths }  = this.props;
    let fileControlProps = {
      file,
      widths,
      callbacks: {removeFile:  this.removeFile.bind(this)}
    };

    return (
      <div className='list-entry-group'>
        <ListEntryTypeColumn { ...{ file, widths } }/>
        <ListEntryFilenameColumn { ...{ file, widths } }/>
        <ListEntryUpdatedColumn { ...{ file, widths } }/>
        <ListEntrySizeColumn { ...{ file, widths } }/>
        <ListEntryControlColumn { ...fileControlProps }/>
      </div>
    );
  }
}

export class ListFolder extends React.Component{
  render(){
    let { file, widths }  = this.props;

    return (
      <div className='list-entry-group'>
        <ListEntryTypeColumn { ...{ is_folder: true, file, widths } }/>
        <ListEntryFilenameColumn { ...{ file, widths } }/>
        <ListEntryUpdatedColumn { ...{ file, widths } }/>
        <ListEntryColumn className='size' widths={widths}/>
        <ListEntryColumn className='control' widths={widths}/>
      </div>
    );
  }
}
