import * as React from 'react';
import { connect } from 'react-redux';
import Icon from './icon';

// this is the Metis control bar, which contains basic operations
// like 'upload file' and 'create folder'

const ControlButton = ({onClick, icon, overlay, title}) => {
  return <button className='control-btn' onClick={ onClick } title={ title }>
    <Icon icon={icon} overlay={overlay}/>
  </button>
}

class ControlBar extends React.Component {
  selectFolder(){
    let { folder_name } = this.props;

    let new_folder_name = prompt("Enter the folder name", "Untitled Folder");
    if (!new_folder_name) return;
    this.props.createFolder(new_folder_name, folder_name);
  }

  selectFile() {
    this.input.click();
  }

  fileSelected(event){
    let { folder_name, fileSelected } = this.props;

    if(event === undefined) return;

    let { files } = this.input;

    for (let i = 0; i < files.length; i++) fileSelected( files[i], folder_name );

    // Reset the input field.
    this.input.value = '';
  }

  render() {
    return (
      <div id='control-bar'>
        <input name='upload-file'
          type='file'
          multiple='multiple'
          style={ {display: 'none'} }
          ref={ (input) => this.input = input }
          onChange={ this.fileSelected.bind(this) }
        />
        <ControlButton onClick={ this.selectFolder.bind(this) } title='Create folder' icon='folder' overlay='plus'/>
        <ControlButton onClick={ this.selectFile.bind(this) } title='Upload file' icon='upload'/>
      </div>
    );
  }
}

export default connect(
  // map state
  null,

  (dispatch) => ({
    fileSelected: (file, folder_name)=>dispatch({ type: 'FILE_SELECTED', file, folder_name }),
    createFolder: (folder_name, parent_folder)=>dispatch({ type: 'CREATE_FOLDER', folder_name, parent_folder })
  })
)(ControlBar);
