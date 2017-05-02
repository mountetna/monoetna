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
    this.spawnWorkers();
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
   * Spawn a worker threads for the uploads.
   */

  spawnWorkers(){

    this['uploadWorker'] = new Worker('./js/workers/uploader.js');
    this['uploadWorker']['onmessage'] = (message)=>{
      
      this.proxyResponse(message);
    };
    this['uploadWorker']['onerror'] = (message)=>{

      console.log(message);
    };

    this['uploadInitializer']= new Worker('./js/workers/upload-initializer.js');
    this['uploadInitializer']['onmessage'] = (message)=>{
      
      this.proxyResponse(message);
    };
    this['uploadInitializer']['onerror'] = (message)=>{

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

        this.initializeFile(action['response']['request']);
        break;
      case 'QUEUE_UPLOAD':

        this.startUpload();
        break;
      case 'PAUSE_UPLOAD':

        var workerMessage = { 'command': 'pause' };
        this['uploadWorker'].postMessage(workerMessage);
        break;
      case 'CANCEL_UPLOAD':

        var delMesg = 'Are you sure you want to remove this upload?'
        if(!confirm(delMesg)) return;
        var workerMessage = { 'command': 'cancel' };
        this['uploadWorker'].postMessage(workerMessage);
        break;
      case 'REMOVE_FILE':

        if(action['fileMetadata'] == undefined) return;
        var delMesg = 'Are you sure you want to remove this file?'
        if(!confirm(delMesg)) return;
        this.removeServerFiles('/remove-file', action['fileMetadata']);
        break;
      case 'REMOVE_FAILED':

        this.removeServerFiles('/remove-failed', action['fileMetadata']);
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
  proxyResponse(message){

    // Set the payload.
    var action = {};
    action['response'] = message['data']['response'];
    action['type'] = message['data']['type'];

    // Dispatch the command.
    this['model']['store'].dispatch(action);

    // Check the queue for another upload or cancel.
    if(message['data']['type'] == 'FILE_UPLOAD_PAUSED') this.startUpload();
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

    var req = ['fileName','originalName','fileSize','projectName','groupName'];
    var reqData = {};

    // Add params not that are needed but not included.
    var state = this['model']['store'].getState();
    reqData['token'] = state['userInfo']['authToken'];

    // Form the upload request.
    req.forEach(function(elem){ reqData[SNAKE_CASE_IT(elem)] = file[elem]; });
    return reqData;
  }

  /*
   * Call to get approval to make an action on Metis.
   */
  requestAuthorization(fileUploadData){

    fileUploadData = this.generateAuthRequest(fileUploadData);
    var request = SERIALIZE_REQUEST(fileUploadData);

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
        'response':response
      };
      this['model']['store'].dispatch(action);
    }
    else{

      console.log('There was an error.');
    }
  }

  initializeFile(authResponse){

    // Get the file from our 'fileUploads' object. 
    var uploadFile = this.getUploadFile(authResponse);
    if(uploadFile == null) return;

    // Normalize the data to send. Add our user token.
    var request = PARSE_REQUEST(uploadFile);
    var state = this['model']['store'].getState();
    request['token'] = state['userInfo']['authToken'];

    // Kick off the 'uploader'.
    var workerMessage = {'command':'init','file':uploadFile,'request':request};
    this['uploadInitializer'].postMessage(workerMessage);
  }

  getUploadFile(authResponse){

    var uploadFile = null;
    var state = this['model']['store'].getState();
    var fileUploads = state['fileData']['fileUploads'];

    fileUploads.forEach(function(elem){

      if((authResponse['groupName'] == elem['groupName']) &&
        (authResponse['projectName'] == elem['projectName']) &&
        (authResponse['fileName'] == elem['fileName'])){

        uploadFile = elem;
      }
    });

    return uploadFile;
  }

  /*
   * For details on the upload/start/pause cycle please refer to the README.md
   * file in the 'workers' folder.
   */
  startUpload(){

    var uploadFile = null;
    var state = this['model']['store'].getState();
    var fileUploads = state['fileData']['fileUploads'];

    for(var a = 0; a < fileUploads['length']; ++a){

      if(fileUploads[a]['status'] == 'queued'){

        uploadFile = fileUploads[a];
      }

      /*
       * If there is an file upload that is active we bail out here and run the
       * pause cycle.
       */
      if(fileUploads[a]['status'] == 'active'){

        var workerMessage = { 'command': 'pause' };
        this['uploadWorker'].postMessage(workerMessage);
        return;
      }
    }

    if(uploadFile == null) return;

    // Normalize the data to send. Add our user token.
    var request = PARSE_REQUEST(uploadFile);
    request['token'] = state['userInfo']['authToken'];

    // Start the upload.
    var workerMessage = {'command':'start','file':uploadFile,'request':request};
    this['uploadWorker'].postMessage(workerMessage);
  }

  removeServerFiles(endPoint, fileMetadata){

    // Serialize the request for POST.
    var request = [];
    for(var key in fileMetadata){
    
      request.push(SNAKE_CASE_IT(key) +'='+ fileMetadata[key]);
    }
    var state = this['model']['store'].getState();
    request.push('token='+ state['userInfo']['authToken']);

    // Request authorization to remove the file.
    AJAX({

      'url': endPoint,
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

      var action = { 'type': 'FILE_REMOVED', 'response': response };
      this['model']['store'].dispatch(action);
    }
    else{

      console.log('There was an error.');
    }
  }

  /*
   * The first step in restarting an upload is to make sure that the last hash, 
   * 'nextBlobHash', that the server saw matches the hash of the file blob on 
   * disk. Once we have make that comparison we also check the file name and
   * size.
   *
   * Next, we calculate the hash of the first 1024 bytes of the file. We send 
   * that info to the server. The server will do the same calculation on it's
   * end. This gives us decent assurance the file on the client matches the
   * partial on the server.
   *
   * Lastly, we recalculate the 'nextBlobHash' of the file on the client.
   * The difference from the first step is that we only hash 1024 bytes. We give
   * this hash to the server. When the server confirms the file upload we reset
   * the 'nextBlobHash' on the server with the new hash. This will reset the
   * upload blob size to 1KiB.
   */

  recoverUpload(uploadFile, fileMetadata){

    var frm = fileMetadata['currentBytePosition'];
    var blob = uploadFile.slice(frm, frm + fileMetadata['nextBlobSize']);
    var fileReader = new FileReader();
    fileReader.onload = function(progressEvent){

      var md5Hash = SparkMD5.ArrayBuffer.hash(this.result);
      metisUploader.checkRecovery(uploadFile, fileMetadata, md5Hash);
    }
    fileReader.readAsArrayBuffer(blob);
  }

  checkRecovery(uploadFile, fileMetadata, md5Hash){

    var fileOk = true;
    if(uploadFile['name'] != fileMetadata['fileName']) fileOk = false;
    if(uploadFile['size'] != fileMetadata['fileSize']) fileOk = false;
    if(fileMetadata['nextBlobHash'] != md5Hash) fileOk = false;

    if(!fileOk){

      alert('The file you selected does not seem to match the record.');
      return;
    }

    this.calcStartingKiB(uploadFile, fileMetadata, md5Hash);
  }

  calcStartingKiB(uploadFile, fileMetadata, md5Hash){

    var blob = uploadFile.slice(0, 1024);
    var fileReader = new FileReader();
    fileReader.onload = function(progressEvent){

      var fkh = SparkMD5.ArrayBuffer.hash(this.result); // firstKiBHash
      metisUploader.calcEndingKiB(uploadFile, fileMetadata, md5Hash, fkh);
    }
    fileReader.readAsArrayBuffer(blob);
  }

  calcEndingKiB(file, fileMetadata, md5Hash, fkh){

    var frm = fileMetadata['currentBytePosition'];
    var blob = file.slice(frm, frm + 1024);
    var fileReader = new FileReader();
    fileReader.onload = function(progressEvent){

      var lkh = SparkMD5.ArrayBuffer.hash(this.result); // lastKiBHash
      metisUploader.recoverUploadOnServer(file,fileMetadata,md5Hash,fkh,lkh);
    }
    fileReader.readAsArrayBuffer(blob);
  }

  recoverUploadOnServer(file, fileMetadata, md5Hash, fristKiBHash, lastKiBHash){

    fileMetadata = this.generateAuthRequest(fileMetadata);
    fileMetadata['first_kib_hash'] = fristKiBHash;
    fileMetadata['last_kib_hash'] = lastKiBHash;
    var request = SERIALIZE_REQUEST(fileMetadata);

    AJAX({

      'url': '/recover-upload',
      'method': 'POST',
      'sendType': 'serial',
      'returnType': 'json',
      'data': request,
      'success': function(response){

        metisUploader.recoverResponse(file, response);
      },
      'error': this['ajaxError'].bind(this)
    });
  }

  recoverResponse(uploadFile, response){

    if(response['success']){

      var action = { 

        'type': 'FILE_UPLOAD_RECOVERED',
        'response': response,
        'uploadFile': uploadFile
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