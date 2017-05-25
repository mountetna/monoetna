import * as React from 'react';

import FileControl from './file-control';

export default class ListEntry extends React.Component{

  constructor(){

    super();
  }

  removeFile(){

    this['props']['callbacks'].removeFile(this['props']['fileInfo']);
  }

  render(){

    var fileInfo = this['props']['fileInfo'];

    var fileControlProps = {

      'fileInfo': fileInfo,
      'callbacks': {'removeFile':  this['removeFile'].bind(this)}
    };

    return (

      <tr className='list-entry-group'>

        <td className='list-entry-icon'>
        </td>
        <td className='list-entry-title-group'>

          <div className='list-entry-file-name'>
            
            {fileInfo['fileName']}
          </div>
          <div className='list-entry-status' title='The current file status.'>

            <span className='light-text'>

              {'uploaded at: '+PARSE_TIMESTAMP(fileInfo['finishTimestamp'])}
            </span>
          </div>
        </td>
        <td className='list-entry-project-group'>

          <div className='list-entry-project-name'>

            {fileInfo['projectName']}
          </div>
          <div className='list-entry-role'>

            <span className='light-text'>

              {fileInfo['role']}
            </span>
          </div>
        </td>
        <td className='list-entry-title-group'>

          <div className='list-entry-file-size'>

            <span className='dark-text' style={{'fontWeight': 900}} >

              {PARSE_BYTES(fileInfo['fileSize'], 1000)}
            </span>
          </div>
          <div className='list-entry-hash'>

            <span className='light-text'>

              {fileInfo['hashingAlgorithm']+': '}
              <span className='mono-text'>
                {fileInfo['hash']}
              </span>
            </span>
          </div>
        </td>
        <FileControl { ...fileControlProps } />
      </tr>
    );
  }
}