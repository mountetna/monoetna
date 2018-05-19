import * as React from 'react';
import * as ReactRedux from 'react-redux';

const ListColumnHead = ({ columnName, widths }) => {
  let columnLabel = columnName.replace(/\b\w/g, l => l.toUpperCase());
  let columnId = `list-${columnName}-column`;

  return <div id={ columnId } className='list-head-title'
    style={ { flexBasis: widths[columnName] } } >
    { columnLabel }
    <div className='list-column-head-arrow-group'>
      <span className='fa fa-chevron-down'></span>
    </div>
  </div>
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
    let { widths } = this.props;
    let fileSelector = {
      id: 'file-selector',
      className: 'file-selector',
      type: 'file',
      name: 'upload-file',
      onChange: this.fileSelected.bind(this)
    };

    return (  
      <div>
        <div id='list-head-group'>
          <ListColumnHead widths={ widths } columnName='type'/>
          <ListColumnHead widths={ widths } columnName='name'/>
          <ListColumnHead widths={ widths } columnName='updated'/>
          <ListColumnHead widths={ widths } columnName='size'/>
          <div id='list-control-column' className='list-head-title'
            style={ { flexBasis: widths.control } }>
            <input { ...fileSelector } />
            <ListControlButton onClick={ this.selectFile.bind(this) } title='Upload file' icon='upload'/>
            <ListControlButton onClick={ this.selectFolder.bind(this) } title='Create folder' icon={ [ 'folder', 'plus' ] }/>
          </div>
        </div>
      </div>
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
