import * as React from 'react';
import { byteFormat, dateFormat } from '../../utils/format';

import FileControl from './file-control';

export default class ListEntry extends React.Component{
  constructor(){
    super();
  }

  removeFile(){
    this.props.callbacks.removeFile(this.props.file);
  }

  render(){
    let { file }  = this.props;
    let fileControlProps = {
      file,
      callbacks: {removeFile:  this.removeFile.bind(this)}
    };

    return (
      <tr className='list-entry-group'>
        <td className='list-entry-icon'/>
        <td className='list-entry-title-group'>
          <div className='list-entry-file-name' title={file.file_name}>
            <a href={file.download_url}>{file.file_name}</a>
          </div>
          <div className='list-entry-status' title='The current file status.'>
            <span className='light-text'>
              {'uploaded at: '+dateFormat(file.finishTimestamp)}
            </span>
          </div>
        </td>
        <td className='list-entry-project-group'>
          <div className='list-entry-project-name'>
            {file.project_name}
          </div>
          <div className='list-entry-role'>
            <span className='light-text'>
              {file.role}
            </span>
          </div>
        </td>
        <td className='list-entry-title-group'>
          <div className='list-entry-file-size'>
            <span className='dark-text' style={{fontWeight: 900}} >
              {byteFormat(file.size, 1000)}
            </span>
          </div>
          <div className='list-entry-hash'>
            <span className='light-text'>
              MD5:
              <span className='mono-text'>
                {file.file_hash}
              </span>
            </span>
          </div>
        </td>
        <FileControl { ...fileControlProps } />
      </tr>
    );
  }
}
