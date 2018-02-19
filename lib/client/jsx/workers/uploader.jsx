/*
 * This class will upload a file as blobs to the Metis upload endpoint.
 */

import md5 from 'md5';
import { postUploadBlob } from '../api/upload_api';
/*
 * In milliseconds, the amount of time to transfer one blob. This ultimately
 * sets the blob size.
 */
const XTR_TIME = 2000;
const MIN_BLOB_SIZE = Math.pow(2, 10); // in bytes
const MAX_BLOB_SIZE = Math.pow(2, 20);

export default (self) => {
  let dispatch = (action) => self.postMessage(action);

  let error = (message) => dispatch(
    { type: 'WORKER_ERROR', worker: 'upload', message }
  );

  // These object will be the data we need to send over the wire.
  let uploadStart = null; // The time at which the upload started in mS.
  let blobWindow = 30; // How many blob uploads to average over.
  let blobUploadTimes = [];

  let pause = false;
  let cancel = false;
  let timeouts = 0; // The number of times an upload has timed out.
  let maxTimeouts = 5; // The number of attepts to upload a blob.

  /*
   * Implementation of a Worker specific function.
   * For message passing into and out of the worker.
   */

  const commandStart = (upload) => {
    let { file, current_byte_position, next_blob_size } = upload;

    if (current_byte_position >= file.size) return;

    let new_blob_size = MAX_BLOB_SIZE;
    let next_byte_position = current_byte_position + next_blob_size;
    let final_byte_position = next_byte_position + new_blob_size;

    // Set the outer limit for 'next_blob_size' and 'endByte'.
    if (final_byte_position > file.size) {
      final_byte_position = file.size;
      new_blob_size = final_byte_position - next_byte_position;
    }

    // Append the current blob to the request.
    let blob = file.slice(current_byte_position, next_byte_position);

    // Hash the next blob.
    let nextBlob = file.slice(next_byte_position, final_byte_position);

    let fileReader = new FileReader();

    fileReader.onload = (event) => {
      let new_blob_hash = md5(fileReader.result);

      let request = { upload, blob, new_blob_size, new_blob_hash };

      sendBlob(request);
    }
    fileReader.readAsArrayBuffer(nextBlob);

    /*
     * This 'response' object usually comes as a response from
     * the server, but on a 'start' command we just fake it to
     * kick off the sequence again.
     */
  }

  let sendBlob = (request) => {
    postUploadBlob(request)
      .then(new_upload => {
        dispatch({ type: 'FILE_UPLOAD_STATUS', upload: new_upload });
        blobComplete(new_upload);
      }).catch(
        () => blobFailed(request)
      );
  }


  self.addEventListener('message', ({data}) => {
    let { upload, command } = data;

    // Check that the incoming data has a valid command.
    switch (command) {
      case 'start':
        commandStart(upload);
        break;
      case 'pause':
        // Set a pause flag.
        commandPause();
        pause = true;
        break;
      case 'cancel':
        // Set a cancel flag.
        cancel = true;
        break;
      default:
        error('Invalid command');
        break;
    }
  });

  /*
   * Calcuates the next blob size based upon upload speed. Also sets the last
   * calculated upload speed on the request.
   */

  let calcNextBlobSize = (current_blob_size) => {
    if (uploadStart != null) {
      let avgTime = averageUploadTime();
      let nextBlobSize = Math.max(
        MIN_BLOB_SIZE,
        Math.min(
          MAX_BLOB_SIZE,
          Math.floor(current_blob_size * (XTR_TIME / avgTime * 0.9))
        )
      );

      /*
       * Make a rough calculation of the upload speed in kilobits per second.
       * This is just for the UI.
       */
      uploadSpeed = (oldBlobSize * 8) / (avgTime / 1000);

      // Reset the timestamp needed for these calculations.
      uploadStart = null;
    }
  };

  let averageUploadTime = () => {
    // Get the current time it took to upload the last blob.
    let uploadTime = (Math.floor(Date.now()) - this.uploadStart);

    // Add the last upload time to an array, limit the array size
    // with a window.
    if (this.blobUploadTimes.length >= this.blobWindow) {
      this.blobUploadTimes.shift();
    }
    this.blobUploadTimes.push(uploadTime);

    // Sum the windowed upload array.
    let uploadTimeSum = this.blobUploadTimes.reduce((a, b)=>{
      return a + b;
    }, 0);

    // Average the upload time over the window.
    return uploadTimeSum / this.blobUploadTimes.length;
  }

  let blobComplete = (response) => {
    switch(response.request.status) {
      case 'active':
        dispatch({ type: 'FILE_UPLOAD_ACTIVE', response });
        /*
         * Check the 'pause' and 'cancel' flags, if set then send
         * the pause or cancel message to the server.
         */
        if(uploader.pause){
          uploader.sendPause(response);
        }
        else if(uploader.cancel){
          uploader.sendCancel(response);
        }
        else{
          uploader.sendBlob(response);
        }
        break;
      case 'paused':
        /*
         * Here the server responed that it got the pause message
         * and the upload is in a paused state. We also clear the
         * 'pause' flag on the uploader.
         */
        uploader.pause = false;
        dispatch({ type: 'FILE_UPLOAD_PAUSED', response });
        break;

      case 'cancelled':
        uploader.cancel = false;
        dispatch({ type: 'FILE_UPLOAD_CANCELLED', response });
        break;

      case 'complete':
        // Send update message back.
        uploader.file = null;
        uploader.request = null;
        dispatch({ type: 'FILE_UPLOAD_COMPLETE', response });
        break;
      default:
        // None
        break;
    }
  }

  let blobFailed = (request) => {
    ++timeouts;
    if(timeouts == maxTimeouts) {
      // Reset the uploader.
      pause = false;
      cancel = false;
      timeouts = 0;

      dispatch({ type: 'FILE_UPLOAD_TIMEOUT', upload });
    }
    else {
      sendBlob(request);
    }
  }
}
