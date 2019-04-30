/*
 * This class will upload a file as blobs to the Metis upload endpoint.
 */

import SparkMD5 from 'spark-md5';
import { postUploadStart, postUploadBlob, postUploadCancel } from '../api/upload_api';
import { setupWorker } from './worker';
const XTR_TIME = 2000; // blob transfer time in ms, determines blob size.
const BLOB_WINDOW = 30; // How many blob uploads to average over.
const MIN_BLOB_SIZE = Math.pow(2, 10); // in bytes
const MAX_BLOB_SIZE = Math.pow(2, 22);
const INITIAL_BLOB_SIZE = Math.pow(2, 10);

const hashBlob = ({nextBlob}) => {
  return new Promise((resolve, reject) => {
    let fileReader = new FileReader();

    fileReader.onload = (event) => {
      let hash = SparkMD5.ArrayBuffer.hash(fileReader.result);

      resolve(hash);
    }

    fileReader.readAsArrayBuffer(nextBlob);
  });
};

const initUpload = (uploader, file, url) => {
  // Hash the next blob for security and error checking.
  uploader.nextBlob = file.slice(0, INITIAL_BLOB_SIZE);

  hashBlob(uploader).then( next_blob_hash => {
    let upload = {
      file_size: file.size,
      next_blob_size: INITIAL_BLOB_SIZE,
      next_blob_hash
    };

    postUploadStart(url, upload)
      .then(uploader.pause)
      .catch(
        () => alert('The upload could not be started.')
      )
  });
};

const cancelUpload = (uploader, upload) => {
  let { status, url, project_name, file_name } = upload;

  if (status == 'complete') {
    uploader.remove(upload);
    return;
  }

  postUploadCancel(url, { project_name, file_name })
    .then(() => uploader.remove(upload))
    .catch(
      (error) => alert('The upload could not be canceled.')
    );
}

// make a new blob based on the current position in the file
const createNextBlob = (uploader, upload) => {
  let { file, current_byte_position, next_blob_size } = upload;

  // after this, we know we have more bytes to send
  if (uploader.paused || current_byte_position >= file.size) return;

  // report the average upload speed, if any
  uploader.setSpeed(upload);

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
        Math.floor(XTR_TIME * uploader.avgSpeed)
      ),
      // but we'll stop at most at here
      MAX_BLOB_SIZE,
      // in fact we should stop when we hit the end of the file
      file.size - next_byte_position
  );

  let final_byte_position = next_byte_position + new_blob_size;

  // Get the two blobs
  let blob_data = uploader.nextBlob;
  uploader.nextBlob = file.slice(next_byte_position, final_byte_position);

  // Hash the next blob.
  hashBlob(uploader).then( new_blob_hash => {
    // Finally, send the request
    let request = {
      blob_data,
      next_blob_size: new_blob_size,
      next_blob_hash: new_blob_hash
    };
    sendBlob(uploader, upload, request);
  });
}

let sendBlob = (uploader, upload, request) => {
  let uploadStart = Date.now();
  let { url } = upload;

  postUploadBlob(url, request)
    .then(new_upload => {
      uploader.addUploadSpeed(uploadStart, upload);
      completeBlob(uploader, {
        ...upload,
        ...new_upload
      });

    }).catch(
      (error) => {
        if (error.fetch) failedBlob(uploader, upload, request);
        else {
          console.log(error);
          throw error
        }
      }
    );
}

const failedBlob = (uploader, upload, request) => {
  if (!uploader.timeout()) sendBlob(uploader, upload, request);
}

const completeBlob = (uploader, new_upload) => {
  let { current_byte_position, file_size } = new_upload;
  if (current_byte_position < file_size) uploader.active(new_upload);
  else uploader.complete(new_upload);
  createNextBlob(uploader, new_upload);
};

export default (self) => {
  let uploader = setupWorker(self, {
    init: ({file,url}) => {
      uploader.paused = false;
      initUpload(uploader, file, url);
    },
    start: ({upload}) => {
      uploader.paused = false;
      createNextBlob(uploader, upload);
    },
    pause: ({upload}) => {
      uploader.paused = true;
      uploader.pause(upload);
    },
    cancel: ({upload}) => {
      uploader.paused = true;
      cancelUpload(uploader, upload);
    }
  });

  uploader.reset = () => {
    uploader.paused = false;
    uploader.cancel = false;
    uploader.timeouts = 0; // The number of times an upload has timed out.
    uploader.maxTimeouts = 5; // The number of attepts to upload a blob.
    uploader.uploadSpeeds = [];
  };

  uploader.reset();

  uploader.status = (upload, status) => uploader.dispatch(
    { type: 'FILE_UPLOAD_STATUS', upload, status }
  );
  uploader.pause = (upload) => uploader.status(upload, 'paused');
  uploader.active = (upload) => uploader.status(upload, 'active');
  uploader.complete = ({file,...upload}) => {
    uploader.status(upload, 'complete');
    uploader.dispatch({ type: 'ADD_FILES', files: [ file ] });
  };
  uploader.remove = (upload) => uploader.dispatch(
    { type: 'FILE_UPLOAD_REMOVED', upload }
  );

  uploader.addUploadSpeed = (uploadStart, upload) => {
    let { next_blob_size } = upload;
    let time = Date.now() - uploadStart;
    uploader.uploadSpeeds.push(next_blob_size / time);
    uploader.uploadSpeeds = uploader.uploadSpeeds.slice(-BLOB_WINDOW);
  }

  // Calcuates the next blob size based upon upload speed.
  uploader.setSpeed = (upload) => {
    let { next_blob_size } = upload;
    let { uploadSpeeds } = uploader;
    if (!uploadSpeeds.length) {
      uploader.avgSpeed = MIN_BLOB_SIZE / XTR_TIME;
      return;
    }

    uploader.avgSpeed = uploadSpeeds.reduce((sum,t)=>sum+t, 0) / uploadSpeeds.length;

    // rough calculation of the upload speed in kbps for the UI
    let upload_speed = uploader.avgSpeed * 8 * 1000;

    uploader.dispatch({ type: 'FILE_UPLOAD_SPEED', upload, upload_speed });
  }

  uploader.timeout = () => {
    ++uploader.timeouts;
    if (uploader.timeouts == uploader.maxTimeouts) {
      uploader.reset();
      uploader.dispatch({ type: 'FILE_UPLOAD_TIMEOUT', upload });
      return false;
    }
    return true;
  }
};
