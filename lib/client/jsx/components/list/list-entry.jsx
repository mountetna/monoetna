import * as React from 'react';
import { byteFormat, dateFormat } from '../../utils/format';

import FileControl from './file-control';
const ListEntryColumn = ({className,widths,children}) =>
  <div className={`list-entry-column-group ${className}`}
    style={ { flexBasis: widths[className] } } >
    { children }
  </div>;

const ListEntryTypeColumn = ({widths}) =>
  <ListEntryColumn className='type' widths={widths}>
    <span className='far fa-file-alt'/>
  </ListEntryColumn>;

const ListEntryControlColumn = ({widths, ...props}) =>
  <ListEntryColumn className='control' widths={widths}>
    <FileControl { ...props } />
  </ListEntryColumn>;

const ListEntrySizeColumn = ({file,widths}) =>
  <ListEntryColumn className='size' widths={widths}>
    <div className='list-entry-file-size'>
      {byteFormat(file.size, 1000)}
    </div>
  </ListEntryColumn>;

const ListEntryUpdatedColumn = ({file, widths}) =>
  <ListEntryColumn className='updated' widths={widths}>
    <div className='list-entry-updated-name'>
      {dateFormat(file.finishTimestamp)}
    </div>
    <div className='list-entry-role'>
      {file.role}
    </div>
  </ListEntryColumn>;

const ListEntryFilenameColumn = ({file, widths}) =>
  <ListEntryColumn className='name' widths={widths}>
    <div className='list-entry-file-name' title={file.file_name}>
      <a href={file.download_url}>{file.file_name}</a>
    </div>
    <div className='list-entry-status' title='The current file status.'>
    </div>
  </ListEntryColumn>;

export default class ListEntry extends React.Component{
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
