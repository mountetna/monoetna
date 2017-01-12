import * as React from 'react'

export default class ListEntry extends React.Component{

  constructor(){

    super();
  }

  render(){

    var fileInfo = this['props']['fileInfo'];
    console.log(fileInfo);

    return (

      <tr className='list-entry-group'>

        <td className='list-entry-icon'>
        </td>
        <td className='list-entry-title-group'>

          <div className='list-entry-file-name'>
            
            { fileInfo['fileName'] }
          </div>
          <div className='list-entry-status' title='The current file status.'>

            <span className='light-text'>

              { 'uploaded at: '+ PARSE_TIMESTAMP(fileInfo['finishTimestamp']) }
            </span>
          </div>
        </td>

        <td className='list-entry-project-group'>

          <div className='list-entry-project-name'>

            { fileInfo['projectName'] }
          </div>
          <div className='list-entry-role'>

            <span className='light-text'>

              { fileInfo['role'] }
            </span>
          </div>
        </td>

        <td className='list-entry-title-group'>

          <div className='list-entry-file-size'>

            <span className='dark-text' style={{ 'fontWeight': 900 }} >

              { PARSE_BYTES(fileInfo['fileSize']) }
            </span>
          </div>
        </td>

        <td className='list-entry-title-group'>

          { /*'controls'*/ }
        </td>

      </tr>
    );
  }
}

/*
      <tr className='list-entry-group'>

        <td className='list-entry-icon'>
        </td>
        <td className='list-entry-title-group'>

          <span className='list-entry-file-name'>
            
            { fileInfo['file_name'] }
          </span>
        </td>
      </tr>
      */