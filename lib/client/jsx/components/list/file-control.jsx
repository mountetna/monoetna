import * as React from 'react';

export default class FileControl extends React.Component{

  constructor(){

    super();
  }

  render(){

    var fileControlBtn = {

      'className': 'file-control-btn',
      'onClick': this['props']['callbacks']['removeFile']
    };

    return (

      <td className='file-control-group'>

        <button { ...fileControlBtn }>

          <span className='glyphicon glyphicon-remove'></span>
        </button>
      </td>
    )
  }
}