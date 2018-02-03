import React from 'react';
import ReactDOM from 'react-dom';
import { Provider } from 'react-redux';
import Cookies from 'js-cookie';

import metisStore from './store';
import MetisUIContainer from './components/metis-ui-container';

import fileData from './reducers/metis-reducer';
import JanusLogReducer from './reducers/janus-log-reducer';
import LastActionReducer from './reducers/last-action-reducer';

import * as authActions from './actions/auth_actions';
import * as fileActions from './actions/file_actions';
import * as uploadActions from './actions/upload_actions';
import * as userActions from './actions/user_actions';

class MetisUploader{
  constructor() {
    this.store = this.createStore();

    this.spawnWorkers();
    this.buildUI();
  }

  createStore() {
    let janusLogReducer = new JanusLogReducer();
    let lastAction = new LastActionReducer();

    // these are (state, action) => new_state
    let reducers = {
      fileData,
      userInfo: janusLogReducer.reducer(),
      lastAction: lastAction.reducer()
    };

    // action handlers to import
    let actions = {
      ...authActions,
      ...fileActions,
      ...uploadActions,
      ...userActions

      // here you may define aliases to other actions,
      // e.g.:
      // returnFile: fileActions.retrieveFile
    };

    return metisStore(reducers, actions);
  }

  initDataStore() {
    // Event hooks from the UI to the Controller
  }

  /*
   * Spawn a worker threads for the uploads.
   */

  spawnWorkers(){
    this.uploadWorker = new Worker('./js/workers/uploader.js');
    this.uploadWorker.onmessage = (message)=>{
      this.proxyResponse(message);
    };
    this.uploadWorker.onerror = (message)=>{
      console.log(message);
    };
    this.uploadInitializer = new Worker('./js/workers/upload-initializer.js');
    this.uploadInitializer.onmessage = (message)=>{
      this.proxyResponse(message);
    };
    this.uploadInitializer.onerror = (message)=>{
      console.log(message);
    };
  }

  buildUI(){
    ReactDOM.render(
      <Provider store={ this.store }>
        <MetisUIContainer />
      </Provider>,
      document.getElementById('ui-group')
    );
  }

  retrieveFiles(){
    let state = this.store.getState();
    let userInfo = state.userInfo;
    if(state.userInfo.loginStatus && !state.userInfo.loginError){
      // Request authorization to upload the file.
      AJAX({
        url: '/retrieve-files',
        method: 'POST',
        sendType: 'serial',
        returnType: 'json',
        data: 'token='+ userInfo.authToken,
        success: this.retrieveFilesResponse.bind(this),
        error: this.ajaxError.bind(this)
      });
    }
  }

  retrieveFilesResponse(response){
    if(response.success){
      let action = { type: 'FILE_METADATA_RECEIVED', fileList: response.file_list };

      this.store.dispatch(action);
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

  /*
   * Responses from the Upload Worker.
   */
  proxyResponse(message){
    // Set the payload.
    let action = { response: message.data.response, type: message.data.type };

    // Dispatch the command.
    this.store.dispatch(action);

    // Check the queue for another upload or cancel.
    if(message.data.type == 'FILE_UPLOAD_PAUSED') this.startUpload();
  }


  generateAuthRequest(file) {
    let req = [
      'fileName',
      'originalName',
      'fileSize',
      'projectName',
      'projectNameFull',
      'groupName'
    ];

    let reqData = {};

    // Add params not that are needed but not included.
    let state = this.store.getState();
    reqData.token = state.userInfo.authToken;

    // Form the upload request.
    req.forEach(elem => { reqData[SNAKE_CASE_IT(elem)] = file[elem]; });
    return reqData;
  }

  initializeFile(authResponse){
    // Get the file from our 'fileUploads' object. 
    let uploadFile = this.getUploadFile(authResponse);
    if(uploadFile == null) return;

    // Normalize the data to send. Add our user token.
    let request = PARSE_REQUEST(uploadFile);
    let state = this.store.getState();
    request.token = state.userInfo.authToken;

    // Kick off the 'uploader'.
    let workerMessage = { command: 'init', file: uploadFile, request };
    this.uploadInitializer.postMessage(workerMessage);
  }

  getUploadFile(authResponse){
    let state = this.store.getState();
    let fileUploads = state.fileData.fileUploads;

    let uploadFile = fileUploads.find(elem => (
      (authResponse.groupName == elem.groupName) &&
      (authResponse.projectName == elem.projectName) &&
      (authResponse.fileName == elem.fileName)
    ))

    return uploadFile;
  }

  /*
   * For details on the upload/start/pause cycle please refer to the README.md
   * file in the 'workers' folder.
   */
  startUpload(){
    let uploadFile = null;
    let state = this.store.getState();
    let { fileUploads } = state.fileData;

    fileUploads.forEach(upload => {
      if (upload.status == 'queued') uploadFile = upload;

      /*
       * If there is an file upload that is active we bail out here and run the
       * pause cycle.
       */
      if(upload.status == 'active'){
        let workerMessage = { command: 'pause' };
        this.uploadWorker.postMessage(workerMessage);
        return;
      }
    })

    if(uploadFile == null) return;

    // Normalize the data to send. Add our user token.
    let request = PARSE_REQUEST(uploadFile);
    request.token = state.userInfo.authToken;

    // Start the upload.
    let workerMessage = {command:'start',file:uploadFile,request};
    this.uploadWorker.postMessage(workerMessage);
  }

  removeServerFiles(endPoint, fileMetadata){
    // Serialize the request for POST.
    let request = [];
    for(let key in fileMetadata){
      request.push(SNAKE_CASE_IT(key) +'='+ fileMetadata[key]);
    }
    let state = this.store.getState();
    request.push('token='+ state.userInfo.authToken);

    // Request authorization to remove the file.
    AJAX({
      url: endPoint,
      method: 'POST',
      sendType: 'serial',
      returnType: 'json',
      data: request.join('&'),
      success: this.removeFileResponse.bind(this),
      error: this.ajaxError.bind(this)
    });
  }

  removeFileResponse(response){
    if(response.success){
      let action = { type: 'FILE_REMOVED', response };
      this.store.dispatch(action);
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
    let frm = fileMetadata.currentBytePosition;
    let blob = uploadFile.slice(frm, frm + fileMetadata.nextBlobSize);
    let fileReader = new FileReader();
    fileReader.onload = function(progressEvent){
      let md5Hash = SparkMD5.ArrayBuffer.hash(this.result);
      metisUploader.checkRecovery(uploadFile, fileMetadata, md5Hash);
    }
    fileReader.readAsArrayBuffer(blob);
  }

  checkRecovery(uploadFile, fileMetadata, md5Hash){
    let fileOk = true;
    if(uploadFile.name != fileMetadata.fileName) fileOk = false;
    if(uploadFile.size != fileMetadata.fileSize) fileOk = false;
    if(fileMetadata.nextBlobHash != md5Hash) fileOk = false;

    if(!fileOk){
      alert('The file you selected does not seem to match the record.');
      return;
    }

    this.calcStartingKiB(uploadFile, fileMetadata, md5Hash);
  }

  calcStartingKiB(uploadFile, fileMetadata, md5Hash){
    let blob = uploadFile.slice(0, 1024);
    let fileReader = new FileReader();
    fileReader.onload = function(progressEvent){
      let fkh = SparkMD5.ArrayBuffer.hash(this.result); // firstKiBHash
      metisUploader.calcEndingKiB(uploadFile, fileMetadata, md5Hash, fkh);
    }
    fileReader.readAsArrayBuffer(blob);
  }

  calcEndingKiB(file, fileMetadata, md5Hash, fkh){
    let frm = fileMetadata.currentBytePosition;
    let blob = file.slice(frm, frm + 1024);
    let fileReader = new FileReader();
    fileReader.onload = function(progressEvent){
      let lkh = SparkMD5.ArrayBuffer.hash(this.result); // lastKiBHash
      metisUploader.recoverUploadOnServer(file,fileMetadata,md5Hash,fkh,lkh);
    }
    fileReader.readAsArrayBuffer(blob);
  }

  recoverUploadOnServer(file, fileMetadata, md5Hash, fristKiBHash, lastKiBHash){
    fileMetadata = this.generateAuthRequest(fileMetadata);
    fileMetadata.first_kib_hash = fristKiBHash;
    fileMetadata.last_kib_hash = lastKiBHash;
    let request = SERIALIZE_REQUEST(fileMetadata);

    AJAX({
      url: '/recover-upload',
      method: 'POST',
      sendType: 'serial',
      returnType: 'json',
      data: request,
      success: (response) => metisUploader.recoverResponse(file, response),
      error: this.ajaxError.bind(this)
    });
  }

  recoverResponse(uploadFile, response){
    if (response.success) {
      let action = { type: 'FILE_UPLOAD_RECOVERED', response, uploadFile };
      this.store.dispatch(action);
    }
    else {
      console.log('There was an error.');
    }
  }

  ajaxError(xhr, config, error){

    console.log(xhr, config, error);
  }
}

// Initilize the class.
let metisUploader = new MetisUploader();
