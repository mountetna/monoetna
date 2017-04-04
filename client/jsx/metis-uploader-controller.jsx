import React from 'react';
import ReactDOM from 'react-dom';
import { Provider } from 'react-redux';

import JanusLogger from './janus-logger-controller';
import MetisModel from './models/metis-model';
import MetisUIContainer from './components/metis-ui-container';

class MetisUploader{

  constructor(){

    this['model'] = null;
    this['uploadWorker'] = null;
    this['janusLogger'] = new JanusLogger();

    this.initDataStore();
    this.spawnUploadWorker();
    this.buildUI();

    /*
     * We pass in the store since we bind events to it. The AJAX callbacks will
     * be dispatched using the store.
     */
    this['janusLogger']['model']['store'] = this['model']['store'];
    this['janusLogger'].checkLog();
  }

  initDataStore(){

    this['model'] = new MetisModel();

    // Event hooks from the UI to the Controller
    this['model']['store'].subscribe(()=>{ 

      var lastAction = this['model']['store'].getState()['lastAction'];
      this.routeAction(lastAction);
    });
  }

  /*
   * Spawn a worker thread for the uploads.
   */
  spawnUploadWorker(){

    this['uploadWorker'] = new Worker('./js/workers/uploader.js');
    this['uploadWorker']['onmessage'] = (message)=>{
      
      this.routeResponse(message);
    };
    
    this['uploadWorker']['onerror'] = (message)=>{

      console.log(message);
    };
  }

  buildUI(){

    ReactDOM.render(

      <Provider store={ this['model']['store'] }>

        <MetisUIContainer />
      </Provider>,
      document.getElementById('ui-group')
    );
  }

  retrieveFiles(){

    var state = this['model']['store'].getState();
    var userInfo = state['userInfo'];
    if(state['userInfo']['loginStatus'] && !state['userInfo']['loginError']){

      // Request authorization to upload the file.
      AJAX({

        'url': '/retrieve-files',
        'method': 'POST',
        'sendType': 'serial',
        'returnType': 'json',
        'data': 'token='+ userInfo['authToken'],
        'success': this['retrieveFilesResponse'].bind(this),
        'error': this['ajaxError'].bind(this)
      });
    }
  }

  retrieveFilesResponse(response){

    if(response['success']){

      var action = { 

        type: 'FILE_METADATA_RECEIVED', 
        fileList: response['file_list'] 
      };

      this['model']['store'].dispatch(action);
    }
    else{

      console.log('There was an error.');
    }
  }

  /*
   * Commands from the UI (via Redux).
   * You would normally see Thunk middleware implemented at the reducer for asyc
   * operations. However, since we are using web workers we don't care about 
   * asyc operations. We receive events from the redux store and we dispatch
   * events to the redux store. This way there is only one entry point to UI, 
   * which is through the redux store, and we do not upset the react/redux 
   * paradigm. Also if there are duplicate actions between here and a reducer.
   * the action in the reducer will run first.
   */
  routeAction(action){

    switch(action['type']){

      case 'AUTHORIZE_FILE':

        if(!this.checkAuthData(action['uploadFile'])){

          alert('The data to upload is not complete.');
          return;
        }
        this.requestAuthorization(action['uploadFile']);
        break;
      case 'FILE_UPLOAD_AUTHORIZED':

        this.initializeFile(action['authResponse']);
        break;
      case 'QUEUE_UPLOAD':

        this.startUpload();
        break;
      case 'PAUSE_UPLOAD':

        this.pauseUpload();
        break;
      case 'CANCEL_UPLOAD':

        var delMesg = 'Are you sure you want to remove this upload?'
        if(!confirm(delMesg)) return;
        this.cancelUpload();
        break;
      case 'REMOVE_FILE':

        if(action['fileMetadata'] == undefined) return;
        var delMesg = 'Are you sure you want to remove this file?'
        if(!confirm(delMesg)) return;
        this.removeFile(action['fileMetadata']);
        break;
      case 'RECOVER_UPLOAD':

        this.recoverUpload(action['uploadFile'], action['fileMetadata']);
        break;
      case 'LOG_IN':

        var email = action['data']['email'];
        var password = action['data']['pass'];
        this['janusLogger'].logIn(email, password);
        break;
      case 'LOGGED_IN':

        this.retrieveFiles();
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
      default:

        //none
        break;
    }
  }

  /*
   * Responses from the Upload Worker.
   */
  routeResponse(message){

    var action = {};
    switch(message['data']['type']){

      case 'initialized':

        action['type'] = 'FILE_INITIALIZED';
        action['initResponse'] = message['data']['response'];
        this['model']['store'].dispatch(action);
        break;
      case 'error':

        //this['uploadWorker'].onerror(message['data']);
        break;
      case 'active':

        action['type'] = 'FILE_UPLOAD_ACTIVE';
        action['uploadResponse'] = message['data']['response'];
        this['model']['store'].dispatch(action);
        break;
      case 'paused':

        action['type'] = 'FILE_UPLOAD_PAUSED';
        action['pauseResponse'] = message['data']['response'];
        this['model']['store'].dispatch(action);

        // Check the queue for another upload or cancel.
        this.startUpload();
        break;
      case 'cancelled':

        /*
         * We could just run the 'removeFile' function from here without 
         * updating the store. However, I am choosing to run an action on the 
         * reducer to update the file metadtata in the store first. Even though 
         * the user will not notice the state change (since it will happen too 
         * fast), this is a more consistent data flow.
         */
        action['type'] = 'FILE_UPLOAD_CANCELLED';
        action['cancelledResponse'] = message['data']['response'];
        this['model']['store'].dispatch(action);

        /*
         * After the store has been updated we can go ahead and extract and 
         * remove the cancelled file(s) from the server.
         */
        var state = this['model']['store'].getState();
        var fileFails = state['fileData']['fileFails'];
        for(var a = 0; a < fileFails['length']; ++a){

          if(fileFails[a]['status'] == 'cancelled'){

            this.removeFile(fileFails[a]);
          }
        }
        
        break;
      case 'complete':

        action['type'] = 'FILE_UPLOAD_COMPLETE';
        action['uploadResponse'] = message['data']['response'];
        this['model']['store'].dispatch(action);
        break;
      default:

        // none
        break;
    }
  }

  checkAuthData(fileUploadData){

    var requestItems = [

      'projectName',
      'projectId',
      'role',
      'fileName',
      'groupId',
      'groupName'
    ];

    var valid = true;
    for(var a = 0; a < requestItems['length']; ++a){

      var key = requestItems[a];
      if(!(key in fileUploadData)) valid = false;  // Checking the key.
      if(fileUploadData[key] == undefined) valid = false; // Checking the value.
    }
    return valid;
  }

  generateAuthRequest(file){

    var fileUploadData = {}

    // Add params not that are needed but not included.
    var state = this['model']['store'].getState();
    fileUploadData['token'] = state['userInfo']['authToken'];

    // Snake case the data for Ruby, isolate the needed params.
    for(var key in file){

      if(key in STATUS_ITEMS){

        var snake_cased_key = SNAKE_CASE_IT(key);
        fileUploadData[snake_cased_key] = file[key];
      }
    }
    return fileUploadData;
  }

  serializeAuthRequset(fileUploadData){

    var request = [];
    for(var key in fileUploadData){

      request.push(key +'='+ fileUploadData[key]);
    }
    return request.join('&');
  }

  /*
   * Call to get approval to make an action on Metis.
   */
  requestAuthorization(fileUploadData){

    fileUploadData = this.generateAuthRequest(fileUploadData);
    var request = this.serializeAuthRequset(fileUploadData);

    // Request authorization to upload the file.
    AJAX({

      'url': '/upload-authorize',
      'method': 'POST',
      'sendType': 'serial',
      'returnType': 'json',
      'data': request,
      'success': this['authorizationResponse'].bind(this),
      'error': this['ajaxError'].bind(this)
    });
  }
  
  /*
   * Route the signed response from Metis or Magma to the web worker.
   */
  authorizationResponse(response){

    if(response['success']){

      var action = {

        'type': 'FILE_UPLOAD_AUTHORIZED',
        'authResponse': response
      };

      this['model']['store'].dispatch(action);
    }
    else{

      console.log('There was an error.');
    }
  }

  initializeFile(authResponse){

    var uploadFile = this.getUploadFile(authResponse);
    if(uploadFile == null) return;

    var initWorker = new Worker('./js/workers/uploader.js');
    initWorker['onmessage'] = (message)=>{
      
      this.routeResponse(message);
    };
    
    initWorker['onerror'] = (message)=>{

      console.log(message);
    };

    var request = PARSE_REQUEST(uploadFile);
    var state = this['model']['store'].getState();
    request['token'] = state['userInfo']['authToken'];

    var workerMessage = {

      'command': 'initialize', 
      'file': uploadFile,
      'request': request
    };

    initWorker.postMessage(workerMessage);
  }

  /*
   * Based upon the dbIndex from the auth response we can extract the actual
   * file object from the redux store.
   */
  getUploadFile(authResponse){

    var request = authResponse['request'];
    var state = this['model']['store'].getState();
    var fileUploads = state['fileData']['fileUploads'];
    var uploadFile = null;
    for(var a = 0; a < fileUploads['length']; ++a){

      if(fileUploads[a]['dbIndex'] == request['dbIndex']){

        uploadFile = fileUploads[a];
        break;
      }
    }

    return uploadFile;
  }

  /*
   * For details on the upload/start/pause cycle please refer to the README.md
   * file in the 'workers' folder.
   */
  startUpload(){

    var state = this['model']['store'].getState();
    var fileUploads = state['fileData']['fileUploads'];
    var uploadFile = null;

    for(var a = 0; a < fileUploads['length']; ++a){

      if(fileUploads[a]['status'] == 'queued'){

        uploadFile = fileUploads[a];
      }

      /*
       * If there is an file upload that is active we bail out here and run the
       * pause cycle.
       */
      if(fileUploads[a]['status'] == 'active'){

        this.pauseUpload();
        return;
      }
    }

    if(uploadFile == null) return;

    var request = PARSE_REQUEST(uploadFile);
    request['token'] = state['userInfo']['authToken'];

    var workerMessage = {

      'command': 'start', 
      'file': uploadFile,
      'request': request
    };

    this['uploadWorker'].postMessage(workerMessage);
  }

  pauseUpload(){

    var workerMessage = { 'command': 'pause' };
    this['uploadWorker'].postMessage(workerMessage);
  }

  cancelUpload(){

    var workerMessage = { 'command': 'cancel' };
    this['uploadWorker'].postMessage(workerMessage);
  }

  recoverUpload(uploadFile, fileMetadata){

    // compare file name
    // compare file size
    // extract blob from upload file
    // hash the blob
    // compare the hash
    console.log(fileMetadata);
  }

  removeFile(fileMetadata){

    // Serialize the request for POST.
    var request = [];
    for(var key in fileMetadata){

      request.push(SNAKE_CASE_IT(key) +'='+ fileMetadata[key]);
    }

    var state = this['model']['store'].getState();
    request.push('authorization_token='+ state['userInfo']['authToken']);
    request.push('signing_algorithm=MD5');

    // Request authorization to remove the file.
    AJAX({

      'url': '/file-remove',
      'method': 'POST',
      'sendType': 'serial',
      'returnType': 'json',
      'data': request.join('&'),
      'success': this['removeFileResponse'].bind(this),
      'error': this['ajaxError'].bind(this)
    });
  }

  removeFileResponse(response){

    if(response['success']){

      var action = {

        'type': 'FILE_REMOVED',
        'oldMetadata': response['old_metadata']
      };

      this['model']['store'].dispatch(action);
    }
    else{

      console.log('There was an error.');
    }
  }

  ajaxError(xhr, config, error){

    console.log(xhr, config, error);
  }
}

// Initilize the class.
var metisUploader = new MetisUploader();