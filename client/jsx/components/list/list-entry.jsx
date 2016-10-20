import * as React from 'react'

export default class ListEntry extends React.Component{

  constructor(){

    super();
  }

  render(){
    
    var fileInfo = this['props']['fileInfo'];
    
    return (

      <div className='list-entry-group'>

        <div className='list-entry-icon'>

          <span className='glyphicon glyphicon-file'></span>
        </div>
        <div className='list-entry-title-group'>
          
          <span className='list-entry-file-name'>
            
            { fileInfo['file_name'] }
          </span>
        </div>
      </div>
    );
  }
}