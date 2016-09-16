/*
 * This class will upload a single file to the Metis upload endpoint.
 */

importScripts('../utils.js');

class SolidUpload{

  /*
   * Implementation of a Worker specific function.
   * For message passing into and out of the worker.
   */
  onmessage(event){

    var message = event.data;

    //Check that the incoming data has a valid command.
    if(!'command' in message){

      postMessage({ message: 'No command.'});
      return;
    }

    var command = message.command;
    switch(command){

      case 'start':

        solidUpload.startUpload(message);
        break;
      case 'stop':

        //Stop the upload.
        break;
      default:

        break;
    }
  }

  /*
   * Take the file input and generates a File Object to send over the wire. The
   * JS File object will set the correct headers and formats for a HTTP POST.
   */
  generateFileObject(file){

    var uploadFile = new FormData();
    uploadFile.append('upload_file', file);
    uploadFile.append('orig_name', file.name);

    return uploadFile;
  }

  /*
  * Pack up the File object with the request data.
  */
  appendRequestToFile(file, request, signature){

    //Unpack the authorized request and reform it for a file upload.
    for(var key in request){

      file.append(key, request[key]);  
    }
    
    //Attach the signature from Magma.
    file.append('signature', signature);

    return file;
  }

  /*
   * This function takes a File oject (which contains the file data) and a 
   * request object that contains the metadata for the upload request.
   */ 
  startUpload(message){

    //Check for the correct data structure.
    if(!'data' in message){

      postMessage({ message: 'No data present.'});
      return;
    }

    if(!'uploadFile' in message.data){

      postMessage({ message: 'No upload file present.'});
      return;
    }

    if(!'request' in message.data){

      postMessage({ message: 'No request data present.'});
      return;
    }

    if(!'signature' in message.data){

      postMessage({ message: 'No signature present.'});
      return;
    }

    var file = message.data.uploadFile;
    var request = message.data.request;
    var signature = message.data.signature;

    //Converts a file input DOM oject to a JS File Object.
    var uploadFile = solidUpload.generateFileObject(file);

    //Appends the request metadata to the File Object for POST upload.
    uploadFile = solidUpload.appendRequestToFile(uploadFile, request, signature);
    
    try{

      AJAX({

        url: '/upload',
        method: 'POST',
        sendType: 'file',
        returnType: 'json',
        data: uploadFile,
        success: this.fileUploadResponse,
        error: this.ajaxError
      });
    }
    catch(error){

      postMessage({ 'type': 'error', message: error.message });
    }
  }

  /*
   * Route the response from Metis.
   */
  fileUploadResponse(response){

    if(response.success){
      
      postMessage({ message: 'File upload complete.' });
    }
    else{

      postMessage({ 'type': 'error', message: 'There was an error.' });
    }
  }

  /*
   * Handle the AJAX errors.
   */
  ajaxError(xhr, ajaxOpt, thrownError){

    postMessage({ xhr: xhr, ajaxOptions: ajaxOpt, thrownError: thrownError });
  }
}

//Initilize the class.
var solidUpload = new SolidUpload();

//Expose the worker specific messaging function.
var onmessage = solidUpload.onmessage;