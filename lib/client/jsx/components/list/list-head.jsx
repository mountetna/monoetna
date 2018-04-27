import * as React from 'react';
import * as ReactRedux from 'react-redux';

class ListHead extends React.Component{
  constructor(){
    super();
  }

  selectFile(){
    /*
     * We are using a button to surragate the file input so we may have 
     * a custom browse button.
     */
    document.getElementById('file-selector').click();
  }

  selectFolder(){
    let folder_name = prompt("Enter the folder name", "Untitled Folder");
    if (!folder_name) return;
    this.props.createFolder(folder_name);
  }

  fileSelected(event){
    if(event === undefined) return;
    let fileSelector = event.target;
    let file = fileSelector.files[0];
    this.props.fileSelected(file);

    // Reset the input field.
    document.getElementById('file-selector').value = '';
  }

  render() {
    let fileSelector = {
      id: 'file-selector',
      className: 'file-selector',
      type: 'file',
      name: 'upload-file',
      onChange: this.fileSelected.bind(this)
    };

    return (  
      <thead>
        <tr id='list-head-group'>
          <th id='list-type-column' className='list-head-title'>
            { 'type ' }
            <div className='list-column-head-arrow-group'>
              <span className='fa fa-caret-down'></span>
            </div>
          </th>
          <th id='list-name-column' className='list-head-title'>
            { 'file name ' }
            <div className='list-column-head-arrow-group'>
              <span className='fa fa-caret-down'></span>
            </div>
          </th>
          <th id='list-project-column' className='list-head-title'>
            { 'project ' }
            <div className='list-column-head-arrow-group'>
              <span className='fa fa-caret-down'></span>
            </div>
          </th>
          <th id='list-size-column' className='list-head-title'>
            { 'size ' }
            <div className='list-column-head-arrow-group'>
              <span className='fa fa-caret-down'></span>
            </div>
          </th>
          <th id='list-control-column' className='list-head-title'>
            <input { ...fileSelector } />
            <button id='file-select-btn' onClick={ this.selectFile.bind(this) }>
              <span className='fa fa-plus white-icon'></span>
              { ' ADD FILE' }
            </button>
            <button id='file-select-btn' onClick={ this.selectFolder.bind(this) }>
              <span className='fa fa-plus white-icon'></span>
              { ' CREATE FOLDER' }
            </button>
          </th>
        </tr>
      </thead>
    );
  }
}

const ListHeadContainer = ReactRedux.connect(
  // map state
  null,

  // map dispatch
  (dispatch) => ({
    fileSelected: (file)=>dispatch({ type: 'FILE_SELECTED', file }),
    createFolder: (folder_name)=>dispatch({ type: 'CREATE_FOLDER', folder_name })
  })
)(ListHead);

export default ListHeadContainer;
