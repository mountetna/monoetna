/*
 * This class will upload a file as blobs to the Metis upload endpoint.
 */

import SparkMD5 from 'spark-md5';
import { postUploadStart, postUploadBlob, postUploadCancel } from '../api/upload_api';
import { setupWorker } from './worker';
import { errorMessage } from '../actions/message_actions';

const XTR_TIME = 2000; // blob transfer time in ms, determines blob size.
const MIN_BLOB_SIZE = Math.pow(2, 10); // in bytes
const MAX_BLOB_SIZE = Math.pow(2, 22);
const INITIAL_BLOB_SIZE = Math.pow(2, 10);

const MAX_UPLOADS = 3; // maximum number of simultaneous uploads

const ZERO_HASH = 'd41d8cd98f00b204e9800998ecf8427e';

const statusFilter = (uploads, s) => uploads.filter(({status}) => status == s);

export default (uploader) => {
  setupWorker(uploader, ({command, ...data}) => {
    if ([ 'start', 'unqueue', 'continue', 'cancel' ].includes(command))
      uploader[command](data);
    else
      console.log('Uploader error', `Invalid command ${command}`);
  });

  Object.assign(uploader, {
    // public commands
    unqueue: ({uploads}) => {
      let active_uploads = statusFilter(uploads, 'active');
      let queued_uploads = statusFilter(uploads, 'queued');

      if (active_uploads > MAX_UPLOADS) {
        // queue some of the active ones
        let overactive_uploads = active_uploads.slice(MAX_UPLOADS);
        overactive_uploads.forEach(uploader.queue);
      } else if (active_uploads < MAX_UPLOADS) {
        // activate some of the queued ones
        let unqueued_uploads = queued_uploads.slice(0, MAX_UPLOADS-active_uploads);
        unqueued_uploads.forEach(uploader.active);
        unqueued_uploads.forEach(upload => uploader.continue({upload}))
      }
    },

    start: ({upload}) => {
      let { file, file_size, file_name, url } = upload;

      let next_blob_size = Math.min(file_size, INITIAL_BLOB_SIZE);

      // Hash the next blob
      let nextBlob = file.slice(0, next_blob_size);

      uploader.hashBlob(nextBlob).then(next_blob_hash => {
        let request = {
          file_size,
          next_blob_size,
          next_blob_hash
        };

        postUploadStart(url, request)
          .then(upload => {
            // this will set the upload status correctly in the upload reducer
            uploader.queue(upload);

            // now we simply broadcast an event for our file:
            uploader.dispatch({ type: 'UNQUEUE_UPLOADS' });
          })
          .catch(
            uploader.warning('Upload failed', error => error)
          )
      });
    },

    continue: ({upload}) => {
      let { file, current_byte_position, next_blob_size, upload_speeds } = upload;

      if (current_byte_position >= file.size) {
        // probably because of a 0 byte file
        let request = {
          action: 'blob',
          blob_data: file.slice(current_byte_position,current_byte_position),
          next_blob_size: 0,
          next_blob_hash: ZERO_HASH,
          current_byte_position
        };
        uploader.sendBlob(upload, request);
        return;
      }

      // after this, we know we have more bytes to send

      // report the average upload speed, if any
      let avgSpeed = !upload_speeds.length
        ? MIN_BLOB_SIZE / XTR_TIME
        : upload_speeds.reduce((sum,t)=>sum+t, 0) / upload_speeds.length;

      // we have already sent up to current_byte_position
      // we will send up to next_byte_position
      let next_byte_position = current_byte_position + next_blob_size;

      // we must also send the hash of the next blob,
      // so we must compute its size, between next_byte_position (the end
      // of this blob) and final_byte_position (the end of the next blob)
      let new_blob_size = Math.min(
          Math.max(
            // we want at least the MIN_BLOB_SIZE
            MIN_BLOB_SIZE,
            // but optimistically, based on the average speed
            Math.floor(XTR_TIME * avgSpeed)
          ),
          // but we'll stop at most at here
          MAX_BLOB_SIZE,
          // in fact we should stop when we hit the end of the file
          file.size - next_byte_position
      );

      let final_byte_position = next_byte_position + new_blob_size;

      // Get the two blobs
      let blob_data = file.slice(current_byte_position, next_byte_position);
      let nextBlob = file.slice(next_byte_position, final_byte_position);

      // Hash the next blob.
      uploader.hashBlob(nextBlob).then( new_blob_hash => {
        // Finally, send the request
        let request = {
          action: 'blob',
          blob_data,
          next_blob_size: new_blob_size,
          next_blob_hash: new_blob_hash,
          current_byte_position
        };

        uploader.sendBlob(upload, request);
      });
    },

    cancel: ({upload}) => {
      let { status, url, project_name, file_name } = upload;

      postUploadCancel(url, { project_name, file_name })
        .then(
          () => uploader.dispatch({ type: 'UPLOAD_FILE_CANCELED', upload })
        )
        .catch(
          uploader.error('Upload cancel failed', error => error)
        );
    },

    // private methods

    reset: () => {
      uploader.timeouts = 0; // The number of times an upload has timed out.
      uploader.maxTimeouts = 5; // The number of attepts to upload a blob.
      uploader.uploadSpeeds = [];
    },

    status: (upload, status) => uploader.dispatch(
      { type: 'UPLOAD_STATUS', upload, status }
    ),
    message: (title, message_handler, message_type) => errorMessage(
      uploader.dispatch, message_type, title, message_handler
    ),
    warning: (title, message) => uploader.message(title, message, 'warning'),
    error: (title, message) => uploader.message(title, message, 'error'),
    notice: (title, message) => uploader.message(title, message, 'notice'),
    pause: (upload) => uploader.status(upload, 'paused'),
    active: (upload) => uploader.status(upload, 'active'),
    complete: (upload) => uploader.status(upload, 'complete'),
    queue: (upload) => uploader.status(upload, 'queued'),

    remove: (upload) => uploader.dispatch({ type: 'REMOVE_UPLOAD', upload }),

    timeout: () => {
      ++uploader.timeouts;
      if (uploader.timeouts == uploader.maxTimeouts) {
        uploader.reset();
        uploader.dispatch({ type: 'UPLOAD_TIMEOUT', upload });
        return false;
      }
      return true;
    },

    hashBlob: (blob) => {
      return new Promise((resolve, reject) => {
        let fileReader = new FileReader();

        fileReader.onload = (event) => {
          let hash = SparkMD5.ArrayBuffer.hash(fileReader.result);

          resolve(hash);
        }

        fileReader.readAsArrayBuffer(blob);
      });
    },

    // make a new blob based on the current position in the file
    sendBlob: (upload, request) => {
      let uploadStart = Date.now();
      let { url } = upload;

      postUploadBlob(url, request)
        .then(new_upload => {
          uploader.dispatch({
            type: 'UPLOAD_SPEED',
            upload_speed: request.blob_data.size / Math.max(1, Date.now() - uploadStart),
            upload
          });
          uploader.completeBlob({
            ...upload,
            ...new_upload
          });
        }).catch(
          error => {
            if (error.fetch) uploader.failedBlob(upload, request);
            else throw error;
          }
        );
    },

    failedBlob: (upload, request) => {
      if (!uploader.timeout()) uploader.sendBlob(upload, request);
    },

    completeBlob: (upload) => {
      let { current_byte_position, file_size } = upload;
      if (current_byte_position < file_size) {
        // update the status
        uploader.status(upload);

        // broadcast that we have finished uploading this blob
        uploader.dispatch({ type: 'UPLOAD_BLOB_COMPLETED', file_name: upload.file_name });
      }
      else {
        // update the status
        uploader.complete(upload);

        // broadcast that we have finished uploading the file
        uploader.dispatch({ type: 'UPLOAD_FILE_COMPLETED', upload });
      }
    }
  });

  uploader.reset();

};
