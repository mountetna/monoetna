/*
 * This class will upload a file as blobs to the Metis upload endpoint.
 */

importScripts('../utils.js');

class BlobUpload{

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

        chunkUpload.extractChunk(message.data.uploadFile);
        break;

      case 'pause':

        //pause the upload
        break;
      case 'stop':

        //Stop the upload.
        break;
      default:

        break;
    }
  }

  extractChunk(file){

    var blob = file.slice(0, 64);
    var uploadBlob = new FormData();
    uploadBlob.append('blob', blob);
    uploadBlob.append('orig_name', file.name);

    try{

      AJAX({

        url: '/upload-blob',
        method: 'POST',
        sendType: 'file',
        returnType: 'json',
        data: uploadBlob,
        success: this.fileUploadResponse,
        error: this.ajaxError
      });
    }
    catch(error){

      postMessage({ 'type': 'error', message: error.message });
    }
  }

  /*
   * Route the response from Metis
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

    // Blob -> ArrayBuffer
    //var uint8ArrayNew  = null;
    //var arrayBufferNew = null;
//    var fileReader     = new FileReader();
//    
//    fileReader.onload  = function(progressEvent) {
//      
//      var uint8Array = new Uint8Array(this.result);
//      console.log(this.result);
//      console.log(uint8Array);

//      var hexString = "";
//      var hexIndex = 0;
//      for(var a = 0; a < uint8Array.byteLength; ++a){

//        var char = uint8Array[a].toString(16).toUpperCase();
//        if(char.length == 1) char = "0"+ char;
//        hexString += " "+char;
//        ++hexIndex;
//        if(hexIndex == 16){

//          console.log(hexString);
//          hexString = "";
//          hexIndex = 0;
//        }
//      }

      
      //arrayBufferNew = this.result;
      //uint8ArrayNew  = new Uint8Array(arrayBufferNew);

      // warn if read values are not the same as the original values
      // arrayEqual from: http://stackoverflow.com/questions/3115982/how-to-check-javascript-array-equals
      //function arrayEqual(a, b) { return !(a<b || b<a); };

      //  if (arrayBufferNew.byteLength !== arrayBuffer.byteLength) // should be 3
      //    
      //    console.warn("ArrayBuffer byteLength does not match");
      //  if (arrayEqual(uint8ArrayNew, uint8Array) !== true) // should be [1,2,3]
      //  
      //    console.warn("Uint8Array does not match");
    //};

    //fileReader.readAsArrayBuffer(blob);

  appendChunkToFile(file, chunk, signature){

  }

  uploadChunk(){

  }
}

//Initilize the class.
var chunkUpload = new ChunkUpload();

//Expose the worker specific messaging function.
var onmessage = chunkUpload.onmessage;