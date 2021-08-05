import * as React from 'react';
import { connect } from 'react-redux';

import ListBody  from './list/list-body';
import ListHead  from './list/list-head';
import FolderBreadcrumb from './folder-breadcrumb';
import ControlBar from './control-bar';

import {selectCurrentFolder} from '../selectors/directory-selector';

const COLUMNS = [
  { name: 'type', width: '60px' },
  { name: 'name', width: '60%' },
  { name: 'status', width: '90px', hide: true },
  { name: 'updated', width: '30%' },
  { name: 'size', width: '10%' },
  { name: 'control', width: '100px', hide: true }
];

const COLUMN_WIDTHS = COLUMNS.reduce( (widths,column) => {
  widths[column.name] = column.width;
  return widths;
}, {} );

const INVALID = '\ninvalid\n';

const InvalidFolder = () => <div className='invalid-folder-view-group'>
  Invalid folder!
</div>;

class FolderView extends React.Component {
  componentDidMount() {
    let { bucket_name, folder_name, retrieveFiles } = this.props;

    retrieveFiles(bucket_name, folder_name);
  }

  selectUpload() {
    this.props.showUploadModal(
      () => this.uploadFileInput.click(),
      () => this.uploadDirInput.click(),
    );
  }

  prepareFiles(event, input) {
    let { bucket_name, folder_name, fileSelected } = this.props;
    if (event === undefined) return;
    let { files } = input;

    for (let i = 0; i < files.length; i++) fileSelected(bucket_name, folder_name, files[i]);

    // Reset the input field.
    input.value = '';
  }

  fileSelected(event){
    this.prepareFiles(event, this.uploadFileInput);
  }

  dirSelected(event) {
    this.prepareFiles(event, this.uploadDirInput);
  }

  selectFolder(){
    let { folder_name, bucket_name, createFolder } = this.props;

    let new_folder_name = prompt("Enter the folder name", "Untitled Folder");

    if (new_folder_name) createFolder(bucket_name, folder_name, new_folder_name);
  }

  selectFolderDownload() {
    let { bucket_name, folder_name } = this.props;
    this.props.listFilesRecursive(bucket_name, folder_name).then(files => {
      return this.props.downloadFilesZip(files, this.props.folder_name, this.props.bucket_name);
    });
  }

  render() {
    let { bucket_name, folder_name, current_folder } = this.props;
    if (current_folder == INVALID) return <InvalidFolder/>;

    let buttons = [
      { onClick: this.selectFolder.bind(this), title: 'Create folder', icon: 'folder', overlay: 'plus', role: 'editor' },
      { onClick: this.selectUpload.bind(this), title: 'Upload file(s)', icon: 'upload', role: 'editor' },
      { onClick: this.selectFolderDownload.bind(this), title: 'Download directory as zip', icon: 'download', role: 'viewer' },
    ];

    return (
      <div className='folder-view-group'>
        <div className='control-group'>
          <FolderBreadcrumb folder_name={folder_name} bucket_name={bucket_name}/>
          <ControlBar buttons={buttons}>
            { /* For uploading individual files */ }
            <input name='upload-file'
              type='file'
              multiple='multiple'
              style={ {display: 'none'} }
              ref={ (input) => this.uploadFileInput = input }
              onChange={ this.fileSelected.bind(this) }
            />

            { /* For uploading directories */ }
            <input name='upload-directory'
              type='file'
              webkitdirectory='webkitdirectory'
              directory='directory'
              multiple='multiple'
              style={ {display: 'none'} }
              ref={ (input) => this.uploadDirInput = input }
              onChange={ this.dirSelected.bind(this) }
            />
          </ControlBar>
        </div>
        <div className='listing-group'>
          <ListHead columns={ COLUMNS }/>
          <ListBody widths={ COLUMN_WIDTHS } folder_name={folder_name} bucket_name={bucket_name}/>
        </div>
      </div>
    );
  }
}

const retrieveFiles = (bucket_name, folder_name) => ({type: 'RETRIEVE_FILES', bucket_name, folder_name});
const fileSelected = (bucket_name, folder_name, file)=>({ type: 'FILE_SELECTED', file, folder_name, bucket_name });
const createFolder = (bucket_name, parent_folder, folder_name)=>({ type: 'CREATE_FOLDER', folder_name, parent_folder, bucket_name });
const showUploadModal = (startFileUpload, startDirectoryUpload) => ({ type: 'SHOW_DIALOG', dialog: { type: 'upload_dialog', startDirectoryUpload, startFileUpload } })
const listFilesRecursive = (bucket_name, folder_name) => ({ type: 'LIST_FILES_RECURSIVE', folder_name, bucket_name });
const downloadFilesZip = (files, folder_name, bucket_name) => ({ type: 'DOWNLOAD_FILES_ZIP', files, folder_name, bucket_name });

export default connect(
  // map state
  state => ({ current_folder: selectCurrentFolder(state) }),

  // map dispatch
  { retrieveFiles, fileSelected, createFolder, showUploadModal, listFilesRecursive, downloadFilesZip }
)(FolderView);
