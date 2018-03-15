import * as React from 'react';

export default class UploadControl extends React.Component{
  constructor(){
    super();
  }

  // Show the start or pause button.
  renderStartPause(){
    let uploadCtrlBtnProps = { className: 'upload-control-btn' };
    let callbacks = this.props.callbacks;

    switch(this.props.upload.status){
      case 'paused':
        uploadCtrlBtnProps.onClick = callbacks.queueUpload;
        return (
          <button { ...uploadCtrlBtnProps }>
            <span className='fa fa-play'></span>
          </button>
        );
      case 'active':
        uploadCtrlBtnProps.onClick = callbacks.pauseUpload;
        return (
          <button { ...uploadCtrlBtnProps }>
            <span className='fa fa-pause'></span>
          </button>
        );
      case 'failed':
        uploadCtrlBtnProps.onClick = callbacks.selectUpload;
        return (
          <button { ...uploadCtrlBtnProps }>
            <span className='fa fa-retweet'></span>
          </button>
        );
        return '';
      default:
        return '';
    }
  }

  render(){
    let uploadControlBtn = {
      className: 'upload-control-btn',
      onClick: this.props.callbacks.cancelUpload
    };

    return (
      <td className='upload-control-group'>
        { this.renderStartPause() }
        <button { ...uploadControlBtn }>
          <span className='fa fa-times'></span>
        </button>
      </td>
    );
  }
}
