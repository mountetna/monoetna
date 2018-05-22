import * as React from 'react';

const UploadButton = ({onClick, icon}) =>
  <button className='upload-control-btn' onClick={onClick}>
    <span className={`fa fa-${icon}`}></span>
  </button>;

export default class UploadControl extends React.Component{
  // Show the start or pause button.
  renderStartPause(){
    let { upload, callbacks } = this.props;

    switch(upload.status){
      case 'paused':
        return <UploadButton
          icon='play'
          onClick={callbacks.queueUpload} />;
      case 'active':
        return <UploadButton
          icon='pause'
          onClick={callbacks.pauseUpload}/>;
      case 'failed':
        return <UploadButton
          icon='retweet'
          onClick={callbacks.selectUpload}/>;
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
