import * as React from 'react';

import PermissionSelector from './permission-selector';
import UploadMeter from './upload-meter';
import UploadControl from './upload-control';

export default class ListUpload extends React.Component{

  constructor(props){

    super(props);

    this['state'] = {

      'componentLock': false,
      'fileNameEditShow': false,
      'fileNameEditActive': false,
      'fileName': this['props']['fileUpload']['fileName']
    };
  }

  componentDidMount(){

    var fileUpload = this['props']['fileUpload'];
    var status = fileUpload['status'];
    if(status == 'authorized' || status == 'active'){

      this.setState({ 

        'componentLock': true ,
        'fileNameEditShow': false,
        'fileNameEditActive': false,
      });
    }
  }

  disabledAlert(){

    alert('You cannot change the file name until the upload is complete.');
  }

  parseFileStatus(){

    var file = this['props']['fileUpload'];
    var status = file['status'];
    var authStamp = file['startTimestamp']
    var date = PARSE_TIMESTAMP(authStamp);
    var user = file['userEmail'];
    switch(file['status']){

      case 'unauthorized':

        status = 'File selected and waiting.'
        break;
      case 'authorized':

        status = 'Auth on '+ date;
        break;
      case 'queued':

        status = 'File upload queued.' 
        break;
      case 'initialized':

        status = 'Initialized by '+ user;
        break;
      case 'active':

        status = 'File uploading...'
        break;
      case 'complete':

        status = 'Uploaded '+ date +' by '+ user;
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

    if(this['state']['fileNameEditActive']){

      return {

        'border': '1px solid #999'
      };
    }
  }

  setInputDisabled(){

    return (this['state']['fileNameEditActive']) ? false : true;
  }

  renderFileNameEditMode(){

    var editBtnProps = {

      'className': 'list-entry-edit-btn',
      'title': 'Edit the file name.',
      'onClick': this['activateFileNameEdit'].bind(this)
    };

    var resetBtnProps = {

      'className': 'list-entry-edit-btn',
      'title': 'Reset the file name.',
      'onClick': this['resetFileName'].bind(this)
    }

    var cancelBtnProps = {

      'className': 'list-entry-edit-btn',
      'title': 'Cancel the file name edit.',
      'onClick': this['deactivateFileNameEdit'].bind(this)
    };

    var saveBtnProps = {

      'className': 'list-entry-edit-btn',
      'title': 'Save the file name.',
      'onClick': this['persistFileName'].bind(this)
    };

    var editShow = this['state']['fileNameEditShow'];
    var editActive = this['state']['fileNameEditActive'];

    if(editShow && !editActive){

      editBtnProps['style'] = { 'display': 'inline-block' };
      cancelBtnProps['style'] = { 'display': 'none' };
      saveBtnProps['style'] = { 'display': 'none' };
      
      var origName = this['props']['fileUpload']['originalName'];
      if(this['state']['fileName'] != origName){

        resetBtnProps['style'] = { 'display': 'inline-block' };
      }
      else{

        resetBtnProps['style'] = { 'display': 'none' };
      }
    }
    else if(editShow || editActive){

      editBtnProps['style'] = { 'display': 'none' };
      resetBtnProps['style'] = { 'display': 'none' };
      cancelBtnProps['style'] = { 'display': 'inline-block' };
      saveBtnProps['style'] = { 'display': 'inline-block' };
    }
    else{

      editBtnProps['style'] = { 'display': 'none' };
      resetBtnProps['style'] = { 'display': 'none' };
      cancelBtnProps['style'] = { 'display': 'none' };
      saveBtnProps['style'] = { 'display': 'none' };
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

    this.setState({ 'fileNameEditShow': true });
  }

  hideFileNameEditMode(event){

    this.setState({ 'fileNameEditShow': false });
  }

  activateFileNameEdit(event){

    if(this['state']['componentLock']){

      this.disabledAlert();
      return;
    }

    this.setState({ 'fileNameEditActive': true });
  }

  deactivateFileNameEdit(event){

    this.setState({ 

      'fileNameEditShow': false,
      'fileNameEditActive': false,
      'fileName': this['props']['fileUpload']['fileName']
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

      var message = 'Not a valid name. Special characters not allowed.\n';
      message += 'Acceptable characters are:  a-z A-Z 0-9 + . - ( ) _';
      alert(message);
      return false;
    }

    return true;
  }

  persistFileName(event){

    if(this['state']['componentLock']){

      this.disabledAlert();
      return;
    }

    // Check for file names.
    
    var newName = this['state']['fileName'];

    if(!this.validateTitle(newName)){

      newName = this['props']['fileUpload']['fileName'];
    }

    this.setState({ 

      'fileNameEditShow': false,
      'fileNameEditActive': false,
      'fileName': newName
    });

    // Bubble the data back to the Redux Store.
    this['props']['fileUpload']['fileName'] = newName;
  }

  updateFileName(event){

    if(this['state']['componentLock']){

      this.disabledAlert();
      return;
    }

    var value = event['target']['value'];
    this.setState({ fileName: value });
  }

  resetFileName(event){

    if(this['state']['componentLock']){

      this.disabledAlert();
      return;
    }

    var origName = this['props']['fileUpload']['originalName'];
    this['props']['fileUpload']['fileName'] = origName;
    this.setState({ fileName: origName });
  }

  persistOnEnter(event){

    event = event || window.event;
    if(event.keyCode == 13 || event.which == 13){

      this.persistFileName(event);
    }
  }

  projectSelected(permission){

    this['props']['fileUpload']['projectName'] = permission['projectName'];
    this['props']['fileUpload']['projectRole'] = permission['role'];
    this['props']['fileUpload']['projectId'] = permission['projectId'];
    this['props']['fileUpload']['groupId'] = permission['groupId'];
  }

  startUpload(){

    this['props']['callbacks'].startFileUpload(this['props']['fileUpload']);
  }

  pauseUpload(){

    //this['props']['callbacks'].pauseUpload();
  }

  cancelUpload(){

    //this['props']['callbacks'].cancelUpload();
  }

  renderPermissionSelector(){

    if(this['props']['fileUpload']['status'] == 'unauthorized'){

      var permissionSelectorProps = {

        'permissions': this['props']['permissions'],
        'fileUpload': this['props']['fileUpload'],
        'editActive': true,
        'callbacks': {

          'projectSelected': this['projectSelected'].bind(this)
        }
      };

      return (

        <PermissionSelector { ...permissionSelectorProps } />
      )
    }
    else{

      var listEntryProjectGroup = {

        'className': 'list-entry-project-group'
      };

      var listEntryStatus = {

        'className': 'list-entry-status list-project-permission-field',
        'title': 'Your project permission for this file.'
      };

      return (

        <td { ...listEntryProjectGroup }>

          <div className='list-project-field'>

            { this['props']['fileUpload']['projectName'] }
          </div>
          <div { ...listEntryStatus }>

            <span className='light-text'>

              { this['props']['fileUpload']['projectRole'] }
            </span>
          </div>
        </td>
      )
    }
  }

  /*
  <td { ...listEntryProjectGroup }>

        <div { ...dropdownGroupProps }>

          <input { ...dropdownInputProps }/>
          <button { ...dropdownBtnProps }>

            <span className='glyphicon glyphicon-triangle-bottom'></span>
          </button>
          <div { ...dropdownTrayProps }>

            { this.addEntries() }
          </div>
        </div>
        <div { ...listEntryStatus }>

          <span className='light-text'>

            { this.setRole() }
          </span>
        </div>
      </td>*/

  render(){

    var listEntryTitleProps = {

      'className': 'list-entry-title-group',
      'onMouseEnter': this['showFileNameEditMode'].bind(this),
      'onMouseLeave': this['hideFileNameEditMode'].bind(this),
    };

    var fileNameInputProps = {

      'className': 'list-entry-file-name list-entry-file-name-input',
      'value': this['state']['fileName'],
      'title': this['state']['fileName'],
      'style': this.setFileNameStyle(),
      'disabled': this.setInputDisabled(),
      'onChange': this['updateFileName'].bind(this),
      'onKeyPress': this['persistOnEnter'].bind(this)
    };

    var uploadControl = {

      'fileUpload': this['props']['fileUpload'],
      'callbacks': {

        'startUpload': this['startUpload'].bind(this),
        'pauseUpload': this['pauseUpload'].bind(this),
        'cancelUpload': this['cancelUpload'].bind(this)
      }
    };

    var listFileStatus = {

      'className': 'list-entry-status',
      'title': 'The current file status.',
      'style': { 'marginTop': '2px' }
    }
    
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
        <UploadMeter fileUpload={ this['props']['fileUpload'] } />
        <UploadControl { ...uploadControl } />
      </tr>
    );
  }
}