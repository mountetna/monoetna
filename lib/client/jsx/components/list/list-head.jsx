import * as React from 'react';
import * as ReactRedux from 'react-redux';

const ListColumnHead = ({ columnName }) => {
  let columnLabel = columnName.replace(/\b\w/g, l => l.toUpperCase());
  let columnId = `list-${columnName}-column`;

  return <th id={ columnId } className='list-head-title'>
    { columnLabel }
    <div className='list-column-head-arrow-group'>
      <span className='fa fa-caret-down'></span>
    </div>
  </th>
};

const ListControlButton = ({onClick, icon, title}) => {
  return <button className='list-control-btn' onClick={ onClick } title={ title }>
    {
      Array.isArray(icon) ?
        <span className='fa-stack fa-fw'>
          <i className={ `fa fa-stack-2x fa-${icon[0]}` }/>
          <i className={ `fa fa-stack-1x fa-${icon[1]} white-icon` }/>
        </span>
      :
        <span className={ `fa fa-fw fa-2x fa-${icon}` }/>
    }
  </button>
}

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
          <ListColumnHead columnName='type'/>
          <ListColumnHead columnName='name'/>
          <ListColumnHead columnName='updated'/>
          <ListColumnHead columnName='size'/>
          <th id='list-control-column' className='list-head-title'>
            <input { ...fileSelector } />
            <ListControlButton onClick={ this.selectFile.bind(this) } title='Upload file' icon='upload'/>
            <ListControlButton onClick={ this.selectFolder.bind(this) } title='Create folder' icon={ [ 'folder', 'plus' ] }/>
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
