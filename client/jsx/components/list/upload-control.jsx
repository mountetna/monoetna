import * as React from 'react'

export default class UploadControl extends React.Component{

  constructor(){

    super();
  }
  
  startUpload(){

    this['props']['callbacks'].startUpload();
  }

  pauseUpload(){

    this['props']['callbacks'].pauseUpload();
  }

  cancelUpload(){

    this['props']['callbacks'].cancelUpload();
  }

  // Show the start or pause button.
  renderStartPause(){

    if(this['props']['fileUpload']['status'] == 'active'){

      return (

        <button className='upload-control-btn' onClick={ this['pauseUpload'].bind(this) }>

          <span className='glyphicon glyphicon-pause'></span>
        </button>
      );
    }
    else{

      return (

        <button className='upload-control-btn' onClick={ this['startUpload'].bind(this) }>

          <span className='glyphicon glyphicon-play'></span>
        </button>
      );
    }
  }

  render(){

    var uploadControlBtn = {

      'className': 'upload-control-btn',
      'onClick': this['cancelUpload'].bind(this)
    };

    return (

      <td className='upload-control-group'>

        { this.renderStartPause() }
        <button { ...uploadControlBtn }>

          <span className='glyphicon glyphicon-remove'></span>
        </button>
      </td>
    );
  }
}