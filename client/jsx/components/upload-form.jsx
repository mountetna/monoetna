import * as React from 'react'

export default class UploadForm extends React.Component{

  constructor(){

    super();
    this.state = {

      textInput: ''
    };
  }

  selectFile(event){

    /**
     * We are using a button to surragate the file input so we may have 
     * a custom browse button.
     */
    document.getElementById('file-selector').click();
  }

  fileSelected(event){

    var fileSelector = event.target;
    var file = fileSelector.files[0];
    this.setState({textInput: file.name});
    this.props.callbacks.fileSelected(file);
  }

  render(){

    return (

      <div id='input-container'>
        
        <input id='text-input' type='text' value={ this.state.textInput }/>
        <input id='file-selector' type='file' name='upload_file' onChange={ (event)=>this.fileSelected(event) }/>
        <button onClick={ (event)=>this.selectFile(event) }>

          browse
        </button>
      </div>
    );
  }
}