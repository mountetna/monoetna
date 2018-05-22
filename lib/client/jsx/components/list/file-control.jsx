import * as React from 'react';

export default class FileControl extends React.Component{
  render(){
    return (
      <div className='file-control-group'>
        <button className='file-control-btn'>
          <span className='glyphicon glyphicon-remove'></span>
        </button>
      </div>
    )
  }
}
