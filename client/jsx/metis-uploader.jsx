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

      type: 'blob',
      user_email: this.userEmail,
      authorization_token: this.authToken,
      original_name: file.name,
      file_size: file.size //in bytes
    };

    this.requestAuthorization(authRequest);
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

      switch(response.request.type){

        case 'solid':

            this.spawnUploadThread(response.request, response.signature);
          break;
        case 'blob':

            this.spwanBlobUploadThread(response.request, response.signature);
          break;
        default:
          break;
      }
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

  /*
   * Spawn a secondary worker thread for the upload.
   */
  spwanBlobUploadThread(request, signature){

    this.blobWorker = new Worker('./js/workers/blob-upload.js');

    this.blobWorker.onmessage = (event)=>{
      
      if(event.data.type == 'error'){
        
        this.blobWorker.onerror(event.data);
      }
      else{

        console.log(event);
      }
    };

    this.blobWorker.onerror = (event)=>{

      console.log(event);
    };

    //Format the upload message for the worker.
    var workerMessage = {

      command: 'start',
      data: {

        uploadFile: this.uploadFile,
        request: request,
        signature: signature
      }
    }

    //Begin the upload.
    this.blobWorker.postMessage(workerMessage);
  }
}

//Initilize the class.
var metisUploader = new MetisUploader();