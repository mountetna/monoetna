import * as React from 'react';

import { byteFormat, dateFormat } from '../../utils/format';
import PermissionSelector from './permission-selector';
import UploadMeter from './upload-meter';
import UploadControl from './upload-control';

export default class ListUpload extends React.Component{

  constructor(props){
    super(props);

    this.state = {
      componentLock: false,
      fileNameEditShow: false,
      fileNameEditActive: false,
      fileName: this.props.upload.file_name
    };
  }

  componentDidMount(){
    let { upload } = this.props;
    let { status } = upload;
    if(status == 'authorized' || status == 'active'){
      this.setState({ 
        componentLock: true ,
        fileNameEditShow: false,
        fileNameEditActive: false,
      });
    }
  }

  disabledAlert(){
    alert('You cannot change the file name until the upload is complete.');
  }

  parseFileStatus(){
    let { upload, user } = this.props;
    let date = dateFormat(upload.startTimestamp);
    let status;

    switch(upload.status){
      case 'paused':
        status = 'waiting for upload.';
        break;

      case 'active':
        status = 'file uploading...';
        break;

      case 'complete':
        status = 'uploaded '+ date +' by '+ user.email;
        break;

      default:
        //none
        break;
    }

    return (
      <span className='light-text'>
        { status }
      </span>
    );
  }

  setFileNameStyle(){
    if(this.state.fileNameEditActive){
      return {
        border: '1px solid #999'
      };
    }
  }

  setInputDisabled(){
    return (this.state.fileNameEditActive) ? false : true;
  }

  renderFileNameEditMode(){
    let editBtnProps = {
      className: 'list-entry-edit-btn',
      title: 'Edit the file name.',
      onClick: this.activateFileNameEdit.bind(this)
    };

    let resetBtnProps = {
      className: 'list-entry-edit-btn',
      title: 'Reset the file name.',
      onClick: this.resetFileName.bind(this)
    }

    let cancelBtnProps = {
      className: 'list-entry-edit-btn',
      title: 'Cancel the file name edit.',
      onClick: this.deactivateFileNameEdit.bind(this)
    };

    let saveBtnProps = {
      className: 'list-entry-edit-btn',
      title: 'Save the file name.',
      onClick: this.persistFileName.bind(this)
    };

    let editShow = this.state.fileNameEditShow;
    let editActive = this.state.fileNameEditActive;

    if(editShow && !editActive){
      editBtnProps.style = { display: 'inline-block' };
      cancelBtnProps.style = { display: 'none' };
      saveBtnProps.style = { display: 'none' };
      
      let origName = this.props.upload.original_name;
      if(this.state.fileName != origName){
        resetBtnProps.style = { display: 'inline-block' };
      }
      else{
        resetBtnProps.style = { display: 'none' };
      }
    }
    else if(editShow || editActive){
      editBtnProps.style = { display: 'none' };
      resetBtnProps.style = { display: 'none' };
      cancelBtnProps.style = { display: 'inline-block' };
      saveBtnProps.style = { display: 'inline-block' };
    }
    else{
      editBtnProps.style = { display: 'none' };
      resetBtnProps.style = { display: 'none' };
      cancelBtnProps.style = { display: 'none' };
      saveBtnProps.style = { display: 'none' };
    }

    return (
      <div className='list-edit-mode-btn-group'>
        <button { ...editBtnProps }>
          <span className='glyphicon glyphicon-pencil'></span>
        </button>
        <button { ...resetBtnProps }>
          <span className='glyphicon glyphicon-refresh'></span>
        </button>
        <button { ...cancelBtnProps }>
          <span className='glyphicon glyphicon-remove'></span>
        </button>
        <button { ...saveBtnProps }>
          <span className='glyphicon glyphicon-ok'></span>
        </button>
      </div>
    );
  }

  showFileNameEditMode(event){
    this.setState({ fileNameEditShow: true });
  }

  hideFileNameEditMode(event){
    this.setState({ fileNameEditShow: false });
  }

  activateFileNameEdit(event){
    if(this.state.componentLock){
      this.disabledAlert();
      return;
    }
    this.setState({ fileNameEditActive: true });
  }

  deactivateFileNameEdit(event){
    this.setState({ 
      fileNameEditShow: false,
      fileNameEditActive: false,
      fileName: this.props.upload.file_name
    });
  }

  validateTitle(newName){
    // Validate that the entry is not blank.
    if(newName == '' || newName == undefined || newName == null){
      alert('Not a valid name. Empty names not allowed.');
      return false;
    }

    // Validate that the entry has no spaces.
    if(/\s/g.test(newName)){
      alert('Not a valid name. Whitespace in names not allowed.');
      return false;
    }

    // Validate that there are no odd characters in the entry.
    if(/[\^\&\'\@\{\}\[\]\,\$\=\!\#\%\+\~]+$/g.test(newName)){
      let message = 'Not a valid name. Special characters not allowed.\n';
      message += 'Acceptable characters are:  a-z A-Z 0-9 + . - ( ) _';
      alert(message);
      return false;
    }

    return true;
  }

  persistFileName(event){
    if(this.state.componentLock){
      this.disabledAlert();
      return;
    }

    // Check for file names.
    
    let newName = this.state.fileName;

    if(!this.validateTitle(newName)){
      newName = this.props.upload.file_name;
    }

    this.setState({ 
      fileNameEditShow: false,
      fileNameEditActive: false,
      fileName: newName
    });

    // Bubble the data back to the Redux Store.
    this.props.upload.fileName = newName;
  }

  updateFileName(event){
    if(this.state.componentLock){
      this.disabledAlert();
      return;
    }

    let value = event.target.value;
    this.setState({ fileName: value });
  }

  resetFileName(event){
    if(this.state.componentLock){
      this.disabledAlert();
      return;
    }

    let origName = this.props.upload.originalName;
    this.props.upload.fileName = origName;
    this.setState({ fileName: origName });
  }

  persistOnEnter(event){
    event = event || window.event;
    if(event.keyCode == 13 || event.which == 13){
      this.persistFileName(event);
    }
  }

  projectSelected(permission) {
    let { upload, callbacks: { selectProject } } = this.props;

    selectProject(upload, permission);
  }

  renderPermissionSelector(){
    let { upload, permissions } = this.props;

    if (upload.status == 'unauthorized') {
      let permissionSelectorProps = {
        permissions,
        upload,
        editActive: true,
        callbacks: {
          projectSelected: this.projectSelected.bind(this)
        }
      };

      return (
        <PermissionSelector { ...permissionSelectorProps } />
      )
    }
    else {
      let listEntryProjectGroup = {
        className: 'list-entry-project-group'
      };

      let permission = permissions.find(p => p.project_name == upload.project_name);

      let listEntryStatus = {
        className: 'list-entry-status list-project-permission-field',
        title: 'Your project permission for this file.'
      };

      return (
        <td { ...listEntryProjectGroup }>
          <div className='list-project-field'>
            { upload.project_name }
          </div>
          <div { ...listEntryStatus }>
            <span className='light-text'>
              { permission.role }
            </span>
          </div>
        </td>
      )
    }
  }

  initializeUpload(){
    this.props.callbacks.initializeUpload(this.props.upload);
  }

  queueUpload(){
    this.props.callbacks.queueUpload(this.props.upload);
  }

  pauseUpload(){
    this.props.callbacks.pauseUpload(this.props.upload);
  }

  cancelUpload(){
    this.props.callbacks.cancelUpload(this.props.upload);
  }

  render(){
    let listEntryTitleProps = {
      className: 'list-entry-title-group',
      onMouseEnter: this.showFileNameEditMode.bind(this),
      onMouseLeave: this.hideFileNameEditMode.bind(this),
    };

    let fileNameInputProps = {
      className: 'list-entry-file-name list-entry-file-name-input',
      value: this.state.fileName,
      title: this.state.fileName,
      style: this.setFileNameStyle(),
      disabled: this.setInputDisabled(),
      onChange: this.updateFileName.bind(this),
      onKeyPress: this.persistOnEnter.bind(this)
    };

    let uploadControl = {
      upload: this.props.upload,
      callbacks: {
        initializeUpload: this.initializeUpload.bind(this),
        queueUpload: this.queueUpload.bind(this),
        pauseUpload: this.pauseUpload.bind(this),
        cancelUpload: this.cancelUpload.bind(this)
      }
    };

    let listFileStatus = {
      className: 'list-entry-status',
      title: 'The current file status.',
      style: { marginTop: '2px' }
    };

    return (
      <tr className='list-entry-group'>
        <td className='list-entry-icon'>
        </td>
        <td { ...listEntryTitleProps }>
          <input { ...fileNameInputProps } />
          { this.renderFileNameEditMode() }
          <div { ...listFileStatus }>
            { this.parseFileStatus() }
          </div>
        </td>
        { this.renderPermissionSelector() }
        <UploadMeter upload={ this.props.upload } />
        <UploadControl { ...uploadControl } />
      </tr>
    );
  }
}
