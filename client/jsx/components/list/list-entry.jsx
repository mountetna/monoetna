import * as React from 'react'

export default class ListEntry extends React.Component{

  constructor(){

    super();
  }

  render(){
    
    var fileInfo = this['props']['fileInfo'];
    
    return (

      <tr className='list-entry-group'>

        <td className='list-entry-icon'>

          <span className='glyphicon glyphicon-file'></span>
        </td>
        <td className='list-entry-title-group'>
          
          <span className='list-entry-file-name'>
            
            { fileInfo['file_name'] }
          </span>
        </td>
      </tr>
    );
  }
}