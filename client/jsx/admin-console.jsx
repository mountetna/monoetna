import React from 'react';
import ReactDOM from 'react-dom';
import { Provider } from 'react-redux';

import JanusLogger from './janus-logger';
import AdminModel from './models/admin-model';
import AdminUIContainer from './components/admin-ui-container';

class AdminConsole{

  constructor(){

    this['model'] = null;
    this['janusLogger'] = new JanusLogger();

    this.initDataStore();

    /*
     * We pass in the store since we bind events to it. The AJAX callbacks will
     * be dispatched using the store.
     */
    this['janusLogger']['model']['store'] = this['model']['store'];
    this['janusLogger'].checkLog();
  }

  initDataStore(){

    this['model'] = new AdminModel();

    // Event hooks from the UI to the Controller
    this['model']['store'].subscribe(()=>{ 

      var lastAction = this['model']['store'].getState()['lastAction'];
      this.routeAction(lastAction);
    });
  }

  buildUI(){

    ReactDOM.render(

      <Provider store={ this['model']['store'] }>

        <AdminUIContainer />
      </Provider>,
      document.getElementById('ui-group')
    );
  }

  routeAction(action){

    switch(action['type']){

      case 'NOT_LOGGED':

        window.location = '/';
        break;
      case 'LOGGED_IN':

        this.buildUI();
        this.fetchAdminData();
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

          'url': METIS_ADDR + endPoint,
          'method': 'POST',
          'sendType': 'serial',
          'returnType': 'json',
          'data': logItems,
          'success': this['adminDataResponse'].bind(this),
          'error': this['ajaxError'].bind(this)
        });
      }
      catch(error){

        //console.log(error);
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

          'url': METIS_ADDR + '/get-permissions',
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

        self.uploadPermissionFile(permissions);
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

    for(var a = 0; a < permissions['length']; ++a){

      if(!('id' in permissions[a]))          permsValid = false;
      if(!('projectId' in permissions[a]))   permsValid = false;
      if(!('projectName' in permissions[a])) permsValid = false;
      if(!('role' in permissions[a]))        permsValid = false;
      if(!('userEmail' in permissions[a]))   permsValid = false;
      if(!('userId' in permissions[a]))      permsValid = false;
    }

    return permsValid;
  }

  uploadPermissionFile(permissions){

    if(COOKIES.hasItem(TOKEN_NAME)){

      //Serialize the request for POST
      var permissionItems = 'token='+ COOKIES.getItem(TOKEN_NAME);

      for(var index in permissions){

        var permission = {};
        for(var key in permissions[index]){

          permission[SNAKE_CASE_IT(key)] = permissions[index][key];
        }

        permissions[index] = permission;
      }

      var encodedPerms = encodeURIComponent(JSON.stringify(permissions));
      permissionItems += '&permissions='+ encodedPerms;

      try{
  
        AJAX({
  
          'url': METIS_ADDR + '/upload-permissions',
          'method': 'POST',
          'sendType': 'serial',
          'returnType': 'json',
          'data': permissionItems,
          'success': this['uploadPermissionResponse'].bind(this),
          'error': this['ajaxError'].bind(this)
        });
      }
      catch(error){
  
        //console.log(error);
      }
    }
  }

  uploadPermissionResponse(response){

    console.log(response);
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
    this.uploadPermissionFile(permission);
/*
    if(COOKIES.hasItem(TOKEN_NAME)){

      //Serialize the request for POST
      var permissionItems = 'token='+ COOKIES.getItem(TOKEN_NAME);
      for(var key in permission){

        permissionItems += '&'+ SNAKE_CASE_IT(key) +'='+ permission[key];
      }

      try{

        AJAX({

          'url': METIS_ADDR + '/save-permission',
          'method': 'POST',
          'sendType': 'serial',
          'returnType': 'json',
          'data': permissionItems,
          'success': this['singlePermissionResponse'].bind(this),
          'error': this['ajaxError'].bind(this)
        });
      }
      catch(error){

        //console.log(error);
      }
    }
*/
  }

//  singlePermissionResponse(response){
//
//    console.log(response);
//  }

  ajaxError(xhr, config, error){

    console.log(xhr, config, error);
  }
}

// Initilize the class.
var adminConsole = new AdminConsole();