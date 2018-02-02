/*
 * This class will initialize an upload to the Metis server from a web client.
 */

importScripts('../libraries/utils.js');
importScripts('../libraries/spark-md5.js');

class UploadInitializer{

  constructor(){

    // These object will be the data we need to send over the wire.
    this.file = null;
    this.request = null;

    this.errorObj=(message)=>({ type:'error', message, response: {} })
    this.timeouts = 0; // The number of times an upload has timed out.
    this.maxTimeouts = 5; // The number of attepts to upload a blob.
  }

  /*
   * Implementation of a Worker specific function.
   * For message passing into and out of the worker.
   */
  onmessage(event){

    var message = event.data;
    uploadInitializer.file = message.file;
    uploadInitializer.request = message.request;

    // Check that the incoming data has a valid command.
    if(!('command' in message)){

      postMessage(uploadInitializer.errorObj('No command.'));
      return;
    }

    if(message.command == 'init'){

      uploadInitializer.initializeUploadSequence();
    }
  }

  /*
   * Create an initial request to the server to initiate an upload. This lets
   * the server know what to expect and also allows the server a chance to
   * resume a previously started upload.
   */
  initializeUploadSequence(){

    /*
     * Here, we are NOT attaching any blob data, but are still using a form
     * for the initialization request. Also keep an eye on that 'SNAKE_CASE_IT'
     * function. The Thin/Ruby server wants it's vars in snake case.
     */
    var initUpReq = new FormData();
    for(var key in  uploadInitializer.request){

      initUpReq.append(SNAKE_CASE_IT(key),  uploadInitializer.request[key]);
    }

    /*
     * Add blob tracking information to the request. This will help us reasemble 
     * and use file upload pausing on the server. All information related to 
     * blobs is in bytes.
     */
    initUpReq.append('current_byte_position', 0);
    initUpReq.append('current_blob_size', 0);
    initUpReq.append('next_blob_size', INITIAL_BLOB_SIZE);

    // Hash the next blob for security and error checking.
    var blob = uploadInitializer.file.slice(0, INITIAL_BLOB_SIZE);
    uploadInitializer.generateBlobHash(blob, initUpReq, this.sendFirstPacket);
  }

  /*
   * Take a blob from the file and hash it.
   */
  generateBlobHash(blob, uploadRequest, callback){

    var fileReader = new FileReader();
    fileReader.onload = function(progressEvent){

      var md5Hash = SparkMD5.ArrayBuffer.hash(this.result);
      uploadRequest.append('next_blob_hash', md5Hash);
      callback(uploadRequest);
    }
    fileReader.readAsArrayBuffer(blob);
  }

  sendFirstPacket(initialUploadRequet){

    try{

      AJAX({

        url: '/upload-start',
        method: 'POST',
        sendType: 'file',
        returnType: 'json',
        data: initialUploadRequet,
        success: uploadInitializer.handleServerResponse,
        error: uploadInitializer.ajaxError
      });
    }
    catch(error){

      postMessage(uploadInitializer.errorObj(error.message));
    }
  }

  handleServerResponse(response){

    uploadInitializer.timeouts = 0; // Reset the timeout counter;

    var errorMessage = 'The server response was malformed.';
    if(!('success' in response) || !('request' in response)){

      postMessage(uploadInitializer.errorObj(errorMessage));
      return;
    }

    response = NORMILIZE_RESPONSE(response);

    if(!response){

      postMessage(uploadInitializer.errorObj(errorMessage));
      return;
    }

    response = TRANSFORM_RESPONSE(response);

    if(!response){

      postMessage(uploadInitializer.errorObj(errorMessage));
      return;
    }

    uploadInitializer.handleResponseRouting(response);  
  }

  handleResponseRouting(response){

    if(response.request.status == 'initialized'){

      postMessage({ type: 'FILE_INITIALIZED', response: response });
    }
  }
}

// Initilize the class.
var uploadInitializer = new UploadInitializer();

// Expose the worker specific messaging function.
var onmessage = uploadInitializer.onmessage;
