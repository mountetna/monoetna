import React from 'react';
import ReactDOM from 'react-dom';
import { Provider } from 'react-redux';

import UserAdminModel from './models/user-admin-model';
import UserAdminViewContainer from './components/user-admin-view-container';
import JanusLoggerController from './janus-logger-controller';

class UserAdminController{

  constructor(){

    this['model'] = null
    this['janusLogger'] = new JanusLoggerController();

    this.initDataStore();

    /*
     * We pass in the store since we bind events to it. The AJAX callbacks will
     * be dispatched using the store.
     */
    this['janusLogger']['model']['store'] = this['model']['store'];
    this['janusLogger'].checkLog();
  }

  initDataStore(){

    this['model'] = new UserAdminModel();

    // Event hooks from the UI to the Controller
    this['model']['store'].subscribe(()=>{ 

      var lastAction = this['model']['store'].getState()['lastAction'];
      this.routeAction(lastAction);
    });
  }

  buildUI(){

    ReactDOM.render(

      <Provider store={ this['model']['store'] }>

        <UserAdminViewContainer />
      </Provider>,
      document.getElementById('ui-group')
    );
  }

  routeAction(action){

    switch(action['type']){

      case 'LOGGED_IN':

        this.buildUI();
        this.fetchAdminData();
        break;
      case 'LOG_OUT':

        this['janusLogger'].logOut(COOKIES.getItem(TOKEN_NAME));
        break;
      case 'LOGGED_OUT':

        window.location = LOGGED_OUT_ADDR();
        break
      case 'NOT_LOGGED':

        window.location = NOT_LOGGED_ADDR();
        break;
      case 'LOGOUT_ALL':

        this.logoutAll();
        break;
      case 'SAVE_PERMISSION':

        this.saveSinglePermission(action['permission']);
        break;
      case 'DOWNLOAD_PERMISSIONS':

        this.downloadPermissions();
        break;
      case 'UPLOAD_PERMISSIONS':

        this.checkPermissionUpload(action['file']);
        break;
      case 'REMOVE_PERMISSION':

        var permissions = [action['permission']];
        this.removePermission(permissions);
        break;
      default:

        // none
        break;
    }
  }

  fetchAdminData(){

    var userInfo = this['model']['store'].getState()['userInfo'];    
    if(!userInfo['masterPerms']){

      window.location = '/';
    }
    else{

      this.adminDataCall('/get-projects');
      this.adminDataCall('/get-permissions');
      this.adminDataCall('/get-users');
    }
  }

  adminDataCall(endPoint){

    if(COOKIES.hasItem(TOKEN_NAME)){

      //Serialize the request for POST
      var logItems = 'token='+ COOKIES.getItem(TOKEN_NAME);

      try{

        AJAX({

          'url': POLYPHEMUS_ADDR + endPoint,
          'method': 'POST',
          'sendType': 'serial',
          'returnType': 'json',
          'data': logItems,
          'success': this['adminDataResponse'].bind(this),
          'error': this['ajaxError'].bind(this)
        });
      }
      catch(error){

        console.log(error);
      }
    }
  }

  adminDataResponse(response){

    console.log(response);
    var action = null;
    if(response['success']){

      if(response.hasOwnProperty('projects')){

        action = { 'type': 'HAS_PROJECTS', 'data': response['projects'] };
      }
      else if(response.hasOwnProperty('permissions')){

        action = { 'type': 'HAS_PERMISSIONS', 'data': response['permissions'] };
      }
      else if(response.hasOwnProperty('users')){

        action = { 'type': 'HAS_USERS', 'data': response['users'] };
      }
      else{

        action = null;
      }
    }

    if(action != null) this['model']['store'].dispatch(action);
  }

  downloadPermissions(){

    if(COOKIES.hasItem(TOKEN_NAME)){

      //Serialize the request for POST
      var logItems = 'token='+ COOKIES.getItem(TOKEN_NAME);

      try{

        AJAX({

          'url': POLYPHEMUS_ADDR + '/get-permissions',
          'method': 'POST',
          'sendType': 'serial',
          'returnType': 'json',
          'data': logItems,
          'success': this['downloadPermsResponse'].bind(this),
          'error': this['ajaxError'].bind(this)
        });
      }
      catch(error){

        //console.log(error);
      }
    }
  }

  downloadPermsResponse(response){

    var object = response['permissions'];
    for(var index in object){

      for(var key in object[index]){

        object[index][CAMEL_CASE_IT(key)] = object[index][key];
        if(key.indexOf('_') != -1) delete object[index][key];
      }
    }
    response['permissions'] = object;

    var json = JSON.stringify(response['permissions'], null, 2);
    var data = 'data:text/json;charset=utf-8,'+ encodeURIComponent(json);
    var downloadAnchor = document.createElement('a');
    downloadAnchor.setAttribute('href', data);
    downloadAnchor.setAttribute('download', 'permissions.json');
    downloadAnchor.click();
  }

  checkPermissionUpload(file){

    var reader = new FileReader();
    var self = this;

    reader.onload = function(progressEvent){

      try{

        var permissions = JSON.parse(this.result);
        if(!self.verifyPermissionUpload(permissions)){

          alert('Your permission file is malformed.');
          return;
        }

        self.uploadPermission(permissions);
      }
      catch(error){

        alert('The file you are attempting to upload is not valid JSON.');
      }
    };

    reader.readAsText(file);
  }

  verifyPermissionUpload(permissions){

    var permsValid = true;

    var permType = Object.prototype.toString.call(permissions);
    if(permType != '[object Array]') permsValid = false;

    var keys = ['projectName','role','userEmail'];
    for(var a = 0; a < permissions['length']; ++a){

      keys.forEach((item)=>{

        if(!(item in permissions[a])) permsValid = false;
      });
    }

    return permsValid;
  }

  uploadPermission(permissions){

    if(COOKIES.hasItem(TOKEN_NAME)){

      //Serialize the request for POST
      var permissionItems = 'token='+ COOKIES.getItem(TOKEN_NAME);
      var encodedPerms = this.normalizePermission(permissions);
      permissionItems += '&permissions='+ encodedPerms;

      try{
  
        AJAX({
  
          'url': POLYPHEMUS_ADDR + '/upload-permissions',
          'method': 'POST',
          'sendType': 'serial',
          'returnType': 'json',
          'data': permissionItems,
          'success': this['permissionsModified'].bind(this),
          'error': this['ajaxError'].bind(this)
        });
      }
      catch(error){
  
        //console.log(error);
      }
    }
  }

  removePermission(permissions){

    if(COOKIES.hasItem(TOKEN_NAME)){

      //Serialize the request for POST
      var permissionItems = 'token='+ COOKIES.getItem(TOKEN_NAME);
      var encodedPerms = this.normalizePermission(permissions);
      permissionItems += '&permissions='+ encodedPerms;

      try{
  
        AJAX({
  
          'url': POLYPHEMUS_ADDR + '/remove-permissions',
          'method': 'POST',
          'sendType': 'serial',
          'returnType': 'json',
          'data': permissionItems,
          'success': this['permissionsModified'].bind(this),
          'error': this['ajaxError'].bind(this)
        });
      }
      catch(error){
  
        //console.log(error);
      }
    }
  }

  permissionsModified(response){

    this.adminDataCall('/get-permissions');
  }

  normalizePermission(permissions){

    for(var index in permissions){

      var permission = {};
      for(var key in permissions[index]){

        permission[SNAKE_CASE_IT(key)] = permissions[index][key];
      }

      permissions[index] = permission;
    }

    return encodeURIComponent(JSON.stringify(permissions));
  }

  /*
   * Here, we make a simple client side check to see if the user and project
   * exisit. We also make this check on the server, but this helps cut down on 
   * unnecessary AJAX calls.
   */
  checkSinglePermission(permission){

    var adminInfo = this['model']['store'].getState()['adminInfo'];

    var users = adminInfo['users'];
    var userExists = false;
    for(var a = 0; a < users['length']; ++a){

      if(users[a]['email'] == permission['userEmail']){

        userExists = true;
      }
    }

    if(!userExists){

      alert("The user specified was not found in the system.");
      return false;
    }

    var projects = adminInfo['projects'];
    var projectExists = false; 
    for(var b = 0; b < projects['length']; ++b){

      if(projects[b]['projectName'].toLowerCase() == permission['projectName']){

        projectExists = true;
      }
    }

    if(!projectExists){

      alert("The project specified was not found in the system.");
      return false;
    }

    return true;
  }

  saveSinglePermission(permission){

    if(!this.checkSinglePermission(permission)) return;

    /*
     * By wrapping the permission in an array we can send it along as a bulk
     * upload...but this 'bulk upload' has only one entry.
     */
    permission = [permission];
    this.uploadPermission(permission);
  }

  logoutAll(){

    var msg1 = 'Are you sure you want to log out all the system users?';
    if(!confirm(msg1)) return;

    var msg2 = 'This action could have a negative impact on the system. \
    Are you really sure?';
    if(!confirm(msg2)) return;

    if(COOKIES.hasItem(TOKEN_NAME)){

      //Serialize the request for POST
      var token = 'token='+ COOKIES.getItem(TOKEN_NAME);
      try{
  
        AJAX({
  
          'url': POLYPHEMUS_ADDR + '/logout-all',
          'method': 'POST',
          'sendType': 'serial',
          'returnType': 'json',
          'data': token,
          'success': function(){ window.location = '/' },
          'error': this['ajaxError'].bind(this)
        });
      }
      catch(error){
  
        window.location = '/';
      }
    }
  }

  ajaxError(xhr, config, error){

    console.log(xhr, config, error);
  }
}

// Initilize the class.
var userAdminController = new UserAdminController();