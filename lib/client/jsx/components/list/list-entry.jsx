import * as React from 'react';
import { byteFormat, dateFormat } from '../../utils/format';

import FileControl from './file-control';

export default class ListEntry extends React.Component{
  constructor(){
    super();
  }

  removeFile(){
    this.props.callbacks.removeFile(this.props.fileInfo);
  }

  render(){
    let fileInfo = this.props.fileInfo;
    let fileControlProps = {
      fileInfo: fileInfo,
      callbacks: {removeFile:  this.removeFile.bind(this)}
    };

    return (
      <tr className='list-entry-group'>
        <td className='list-entry-icon'/>
        <td className='list-entry-title-group'>
          <div className='list-entry-file-name' title={fileInfo.fileName}>
            {fileInfo.fileName}
          </div>
          <div className='list-entry-status' title='The current file status.'>
            <span className='light-text'>
              {'uploaded at: '+dateFormat(fileInfo.finishTimestamp)}
            </span>
          </div>
        </td>
        <td className='list-entry-project-group'>
          <div className='list-entry-project-name'>
            {fileInfo.projectName}
          </div>
          <div className='list-entry-role'>
            <span className='light-text'>
              {fileInfo.role}
            </span>
          </div>
        </td>
        <td className='list-entry-title-group'>
          <div className='list-entry-file-size'>
            <span className='dark-text' style={{fontWeight: 900}} >
              {byteFormat(fileInfo.fileSize, 1000)}
            </span>
          </div>
          <div className='list-entry-hash'>
            <span className='light-text'>
              {fileInfo.hashingAlgorithm+': '}
              <span className='mono-text'>
                {fileInfo.hash}
              </span>
            </span>
          </div>
        </td>
        <FileControl { ...fileControlProps } />
      </tr>
    );
  }
}
