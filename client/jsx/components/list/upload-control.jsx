import * as React from 'react'

export default class UploadControl extends React.Component{

  constructor(){

    super();
  }
  
  startUpload(){

    //this['props'].startUpload(this['props']['fileUpload']);
  }

  pauseUpload(){

    //this['props'].pauseUpload(this['props']['fileUpload']);
  }

  cancelUpload(){

    //this['props'].cancelUpload(this['props']['fileUpload']);
  }

  // Show the start or pause button.
  renderStartPause(){

    var fileUpload = this['props']['fileUpload'];
    
    if(fileUpload['status'] == 'active'){

      return (

        <button className='upload-control-btn' onClick={ this.pauseUpload.bind(this) }>

          <span className='glyphicon glyphicon-pause'></span>
        </button>
      );
    }
    else{

      return (

        <button className='upload-control-btn' onClick={ this.startUpload.bind(this) }>
          
          <span className='glyphicon glyphicon-play'></span>
        </button>
      );
    }
  }

  render(){

    return (

      <td className='upload-control-group'>

        { this.renderStartPause() }
        <button className='upload-control-btn' onClick={ this.cancelUpload.bind(this) }>

          <span className='glyphicon glyphicon-remove'></span>
        </button>
      </td>
    );
  }
}