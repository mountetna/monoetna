import React from 'react';
import ReactDOM from 'react-dom';

import UploadForm from './components/upload-form';

class MetisUploader{

  constructor(){

    this.buildUI();
    this.uploadFile = null;
    this.uploadWorker = null;

    //these two items are dummy entries
    this.userEmail = 'jason.cater@ucsf.edu';
    this.authToken = 'blahf8cfc63531b6ed753bd536f6e12d578c';
  }

  buildUI(){

    //an object with functions so the component can call back to the controller
    var callbacks = {};
    callbacks.fileSelected = this.fileSelected.bind(this);

    ReactDOM.render(

      <UploadForm callbacks={ callbacks } />,
      document.getElementById('ui-container')
    );
  }

  /*
   * File selected, now call Magma to get an HMAC Auth Message
   */
  fileSelected(file){

    //save the file handle
    this.uploadFile = file;

    var authRequest = {

      type: 'upload',
      auth_token: this.authToken,
      user_email: this.userEmail,
      file_name: file.name
    };

    this.requestAuthorization(authRequest);
    //this.spawnChunkUploadThread();
  }

  /*
   * Call Magma to get approval to make an action on Metis
   * type = PUT, GET, DEL
   */
  requestAuthorization(authRequest){

    //Serialize the request for POST
    var request = [];
    for(var key in authRequest){

      request.push(key +'='+ authRequest[key]);
    }
    request = request.join('&');

    AJAX({

      url: './magma-end-point', //replace this URL with the Magma End Point
      method: 'POST',
      sendType: 'serial',
      returnType: 'json',
      data: request,
      success: this.authorizationResponse.bind(this),
      error: this.ajaxError
    });
  }
  
  /*
   * Route the signed response from Magma to the web worker.
   */
  authorizationResponse(response){

    if(response.success){

      this.spawnUploadThread(response.request, response.signature);
    }
    else{

      console.log('There was an error.');
    }
  }

  /*
   * Spawn a secondary worker thread for the upload.
   */
  spawnUploadThread(request, signature){

    //Setup the worker.
    this.uploadWorker = new Worker('./js/workers/solid-upload.js');

    this.uploadWorker.onmessage = (event)=>{
      
      if(event.data.type == 'error'){

        this.uploadWorker.onerror(event.data);
      }
      else{

        console.log(event.data.message);
      }
    };

    this.uploadWorker.onerror = (event)=>{

      if(event.type === 'error') console.log(event.message);      
    };

    //Format the upload message for the worker.
    var workerMessage = {

      command: 'start',
      data: {

        uploadFile: this.uploadFile,
        request: request,
        signature: signature
      }
    };

    //Begin the upload.
    this.uploadWorker.postMessage(workerMessage);
  }

  spawnChunkUploadThread(){

    this.chunkWorker = new Worker('./js/workers/chunk-upload.js');

    this.chunkWorker.onmessage = (event)=>{
      
      console.log(event);
    };

    this.chunkWorker.onerror = (event)=>{

      console.log(event);
    };

    //Format the upload message for the worker.
    var workerMessage = {

      command: 'start',
      data: {

        uploadFile: this.uploadFile,
        //signature: signature
      }
    }

    this.chunkWorker.postMessage(workerMessage);
  }
}

//Initilize the class.
var metisUploader = new MetisUploader();