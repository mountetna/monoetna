import * as React from 'react';

export default class UploadControl extends React.Component{

  constructor(){

    super();
  }

  // Show the start or pause button.
  renderStartPause(){

    var uploadCtrlBtnProps = { 'className': 'upload-control-btn' };
    var callbacks = this['props']['callbacks'];

    switch(this['props']['fileUpload']['status']){

      case 'unauthorized':

        uploadCtrlBtnProps['onClick'] = callbacks['initializeUpload'];
        return (

          <button { ...uploadCtrlBtnProps }>

            <span className='glyphicon glyphicon-arrow-right'></span>
          </button>
        );
      case 'initialized':
      case 'paused':

        uploadCtrlBtnProps['onClick'] = callbacks['queueUpload'];
        return (

          <button { ...uploadCtrlBtnProps }>

            <span className='glyphicon glyphicon-play'></span>
          </button>
        );
      case 'active':

        uploadCtrlBtnProps['onClick'] = callbacks['pauseUpload'];
        return (

          <button { ...uploadCtrlBtnProps }>

            <span className='glyphicon glyphicon-pause'></span>
          </button>
        );
      case 'failed':

        uploadCtrlBtnProps['onClick'] = callbacks['selectUpload'];
        return (

          <button { ...uploadCtrlBtnProps }>

            <span className='glyphicon glyphicon-retweet'></span>
          </button>
        );
        return '';
      default:

        return '';
    }
  }

  render(){

    var uploadControlBtn = {

      'className': 'upload-control-btn',
      'onClick': this['props']['callbacks']['cancelUpload']
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