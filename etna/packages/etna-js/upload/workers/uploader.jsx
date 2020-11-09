import 'regenerator-runtime/runtime';
import SparkMD5 from 'spark-md5';
import {postUploadStart, postUploadBlob, postUploadCancel} from '../api/upload_api';
import {setupWorker} from './worker';
import {checkStatus} from "../../utils/fetch";

import {fileKey} from "../../utils/file";
import {Cancellable} from "../../utils/cancellable";
import {Debouncer} from "../../utils/debouncer";
import {Schedule} from "../../utils/schedule";
import {assertIsSome} from "../../utils/asserts";
import {Subscription} from "../../utils/subscription";

const XTR_TIME = 2000; // blob transfer time in ms, determines blob size.
const MIN_BLOB_SIZE = Math.pow(2, 10); // in bytes
const MAX_BLOB_SIZE = Math.pow(2, 22);
const MAX_UPLOADS = 3; // maximum number of simultaneous uploads
const UPLOAD_SPEED_WINDOW = 30;
export const ZERO_HASH = 'd41d8cd98f00b204e9800998ecf8427e';

// Commands in
export const ADD_UPLOAD = 'ADD_UPLOAD';
export const CANCEL_UPLOAD = 'CANCEL_UPLOAD';
export const PAUSE_UPLOAD = 'PAUSE_UPLOAD';

// Events out
export const UPDATE_UPLOADS = 'UPDATE_UPLOADS';
export const UPLOAD_COMPLETE = 'UPLOAD_COMPLETE';
export const UPLOAD_ERROR = 'UPLOAD_ERROR';

// Main entrypoint for the worker.  Creates uploader and connects it to the postMessage and message event handlers to
// communicate with host webpage.
export default (worker, opts = {}) => {
  const {dispatch, subscribe} = setupWorker(worker);
  const uploader = new Uploader({dispatch}, {subscribe}, opts);
  uploader.start();
  return uploader;
};

// Event posted by the uploader when any upload state changes, includes bundle of all upload objects.
export function UpdateUploadsEvent(uploads) {
  return {
    type: UPDATE_UPLOADS,
    uploads: uploads,
  };
}

// Event posted by the uploader when any upload successfully completes.
export function UploadCompleteEvent(upload) {
  return {
    type: UPLOAD_COMPLETE,
    upload,
  }
}

// Event posted by the uploader for any user facing failures during upload.
export function UploadErrorEvent(message_type, title, message, upload = undefined) {
  return {
    type: UPLOAD_ERROR,
    message: message,
    title: title,
    message_type: message_type,
    upload: upload,
  };
}

// Command to send to the uploader that adds a new / restarts an existing upload
// by the given fileKey of that upload object.
export function AddUploadCommand(upload) {
  return {
    command: ADD_UPLOAD,
    upload: upload,
  };
}

// Command to send to the uploader that cancels and removes an existing upload
// by the given fileKey of that upload object.  Safe if the upload does not exist.
export function CancelUploadCommand(upload) {
  return {
    command: CANCEL_UPLOAD,
    upload: upload,
  };
}

// Command to send the uploader that pauses an existing upload by the given fileKey of
// that upload object.  Has no effect if the given upload did not already exist.
export function PauseUploadCommand(upload) {
  return {
    command: PAUSE_UPLOAD,
    upload: upload,
  };
}

export function Upload({
    file_name,
    project_name,
    url,
    file,
    next_blob_size = 0,
    upload_speeds = [],
    status = 'preparing',
    current_byte_position = 0,
    next_blob_hash = ZERO_HASH,
    file_size = file.size,
    metis_uid,
    attribute_name,
    model_name,
    record_name,
    original_filename
  }) {
    assertIsSome({ file_name, project_name, url, file });

  return {
    status: status,
    file_name: file_name,
    project_name: project_name,
    next_blob_size: next_blob_size,
    next_blob_hash: next_blob_hash,
    upload_speeds: upload_speeds,
    current_byte_position: current_byte_position,
    file: file,
    url: url,
    file_size: file_size,
    metis_uid,
    attribute_name,
    model_name,
    record_name,
    original_filename
  }
}

// Just convenience for the fact that 'unpausing' is really identical to adding a paused upload.
export const UnpauseUploadCommand = AddUploadCommand;

// Exported for test purposes
export function hashBlob(blob) {
  return new Promise((resolve, reject) => {
    let fileReader = new FileReader();

    fileReader.onload = (event) => {
      let hash = SparkMD5.ArrayBuffer.hash(fileReader.result);

      resolve(hash);
    }

    fileReader.readAsArrayBuffer(blob);
  });
}

// Main Service Class for managing uploads.  Uploader maintains its own uploads state
// and performs side effects to push the data to the server.  The main interface is through
// the eventSink and eventSource objects passed to the constructor.
export class Uploader {
  constructor(
    // receives events from the uploader.  See Event classes above.
    eventSink = {
      dispatch(data) {
      }
    },
    eventSource = {
      // main input interface.  After Uploader.start is called, it passes a handler to this method that is
      // expected to provide that handler future commands.  See Command classes above.
      subscribe(f) {
      }
    },
    {
      maxUploads = MAX_UPLOADS,
      minBlobSize = MIN_BLOB_SIZE,
      maxBlobSize = MAX_BLOB_SIZE,
      maxUploadFailures = 5,
      backoffFactor = 1000,
      debouncerOptions = {},
    } = {},
  ) {
    this.numUploadFailuresByKey = {};
    this.maxUploadFailures = maxUploadFailures;
    this.backOffFactor = backoffFactor;

    this.eventSource = eventSource;
    this.eventSink = eventSink;

    // Maps uploadKeys -> upload object state
    this.uploads = {};
    this.maxUploads = maxUploads;
    this.minBlobSize = minBlobSize;
    this.maxBlobSize = maxBlobSize;

    // Maps uploadKeys -> cancellables
    this.running = {};

    // Keeps track of pending synchronize work, useful for tests to await all work being complete.
    this.schedule = new Schedule();

    // Batches updates together into groups that do not occur more than once every 100ms.
    this.debouncer = new Debouncer({maxGating: 100, windowMs: 100, eager: true, ...debouncerOptions});

    this.subscription = new Subscription();
  }

  start() {
    // Idempotent start.
    this.start = () => null;

    this.eventSource.subscribe(({command, upload}) => {
      switch (command) {
        case ADD_UPLOAD:
          this.addUpload(AddUploadCommand(upload));
          break;
        case CANCEL_UPLOAD:
          this.cancelUpload(CancelUploadCommand(upload));
          break;
        case PAUSE_UPLOAD:
          this.pauseUpload(PauseUploadCommand(upload));
          break;
        default:
          console.error('Unexpected uploader command', command);
      }
    });
  }

  cancelRunningUpload(key) {
    let running = this.running[key];
    delete this.running[key];
    running.cancel();
  }

  restartUploadCancellable(key) {
    if (key in this.running) this.cancelRunningUpload(key);
    const newCancellable = this.running[key] = new Cancellable();
    this.subscription.addCleanup(newCancellable.cancel);
    return newCancellable;
  }

  // Called anytime the uploads state is changed.
  // Cancels any uploads that are no longer in the active state,
  // then attempts to find an upload to start if there is less than maxUploads in active or paused state.
  // Recursive -- when initiating a new upload, this method will re-enter via updateUpload, thus
  // each turn about will activate one upload at a time until no new ones need activation.
  synchronize() {
    // Cancel any running that are no longer in uploads or are paused
    for (let key in this.running) {
      // Cancel any non active uploads that are still running.  This immediately stops
      // those running to end via the cancellable.
      if (!(key in this.uploads) || this.uploads[key].status !== 'active') {
        this.cancelRunningUpload(key);
      }
    }

    // Count the number of uploads that are either active or paused; they count against the maxUploads count
    const runningOrPaused = Object.values(this.uploads).filter(({status}) => status === 'active' || status === 'paused');
    let availableStartCount = this.maxUploads - runningOrPaused.length;

    // Start any queued uploads up to the max
    for (let key in this.uploads) {
      // Do not continue if there are not more available start slots.
      if (availableStartCount <= 0) break;
      const upload = this.uploads[key];
      if (upload.status === 'queued') {
        this.updateUpload({...this.uploads[key], status: 'active'}, false);
        this.startUpload(key);
        availableStartCount -= 1;
      }
    }

    this.schedule.addWork(
      this.debouncer.ready(this.subscription.addSubscribedCallback(() =>
        this.eventSink.dispatch(UpdateUploadsEvent(this.uploads))))
    );
  }

  startUpload(uploadKey, reset = false) {
    let {file_name, status} = this.uploads[uploadKey];
    // Invariant control, to ensure that running uploads have already been marked as active.
    if (status !== 'active') {
      throw new Error(`Cannot restart download of ${file_name}, was ${status} and not active`);
    }

    const cancellable = this.restartUploadCancellable(uploadKey);
    return this.schedule.addWork(cancellable.run(this.uploadLoop(uploadKey, reset))
      .then(
        ({cancelled, result: upload}) => {
          if (cancelled) return {cancelled}
          upload = this.updateUpload({...upload, status: 'complete'});
          this.eventSink.dispatch(UploadCompleteEvent(upload));
        },
        err => this.reportErrorAndPause(uploadKey, `Upload of file ${file_name}`, err))).catch(e => console.warn(e));
  }

  * uploadLoop(uploadKey, reset) {
    let {file, file_size, url} = this.uploads[uploadKey];
    let next_blob_size = Math.min(file_size, this.minBlobSize);
    let nextBlob = file.slice(0, next_blob_size);

    const next_blob_hash = yield hashBlob(nextBlob);

    let request = {
      file_size,
      next_blob_size,
      next_blob_hash,
      reset
    };

    let {current_byte_position} = this.updateUpload(yield postUploadStart(url, request));
    let unsentZeroByteFile = current_byte_position === 0;

    // Continue uploading while either there is more to send, or nothing has been sent (case of 0 byte file)
    while (current_byte_position < file.size || unsentZeroByteFile) {
      ({current_byte_position} = yield* this.sendNextBlob(uploadKey));

      unsentZeroByteFile = false;
    }

    return this.uploads[uploadKey];
  }

  * sendNextBlob(uploadKey) {
    const upload = this.uploads[uploadKey];
    const {url, upload_speeds} = upload;

    const request = yield this.prepareNextBlobParams(upload);
    const uploadStart = Date.now();
    const {response, err} = yield postUploadBlob(url, request, false).then(response => ({response}), err => ({err}));

    const serverUpload = yield* this.checkForRetry(uploadKey, response, err);

    const uploadSpeed = request.blob_data.size / Math.max(1, Date.now() - uploadStart);

    return this.updateUpload({
      ...serverUpload,
      upload_speeds: [...upload_speeds, uploadSpeed].slice(-UPLOAD_SPEED_WINDOW)
    });
  }

  * checkForRetry(uploadKey, response = null, err = null) {
    let reset = false;
    let retry = false;

    // Network errors won't produce a response
    if (response == null) {
      retry = true;
    } else if (response.status === 422) {
      // The blob has a mismatch on its hash or content.  It's possible that a concurrent upload may have occurred
      // and corrupted the server's view of the data.
      reset = true;
      retry = true;
    } else if (response.status >= 500) {
      retry = true;
    }

    // Retries will not be attempted, however, if the number of retries has already succeeded the max.
    if (retry && ++this.numUploadFailuresByKey[uploadKey] >= this.maxUploadFailures) {
      retry = false;
    }

    if (!retry) {
      // Handle the case that the deserialization of the response fails during body transfer or
      // due to non json response.
      if (response) {
        const bodyPromise = checkStatus(response);
        const {body, err} = yield bodyPromise.then(body => ({body}), err => ({err}));
        if (err) {
          return yield* this.checkForRetry(uploadKey, null, err);
        }

        return body;
      }
      else throw err;
    }

    let retryCount = this.numUploadFailuresByKey[uploadKey];
    console.warn('Retrying request for', uploadKey, 'failed with', response ? response.status : err)

    const timeout = Math.pow(2, retryCount) * this.backOffFactor;

    yield new Promise((resolve, reject) => {
      // 2000ms
      // 4000ms
      // 8000ms
      // 16000ms
      setTimeout(() => {
        // this will cancel the current cancellable and start a new processing of the upload.
        try {
          this.startUpload(uploadKey, reset);
          resolve();
        } catch(e) {
          console.error(e);
          reject(e);
        }
      }, timeout);
    });
  }

  reportErrorAndPause(uploadKey, title, err) {
    let upload = this.uploads[uploadKey];
    upload = this.pauseUpload({upload});
    return this.reportError(title, err, upload);
  }

  reportError(title, err, upload = undefined) {
    const self = this;

    if (err instanceof Promise) {
      return err.then(
        ({error}) => reportErrorAndThrow({title, message: error, e: new Error(error)}),
        e => reportErrorAndThrow({title, e}),
      );
    }

    reportErrorAndThrow({title});

    function reportErrorAndThrow({title, message = 'An unknown error occurred, try again later', e = err}) {
      self.eventSink.dispatch(UploadErrorEvent('error', title, message, upload));
      throw e;
    }
  }

  prepareNextBlobParams(upload) {
    const {file, current_byte_position, next_blob_size, upload_speeds, metis_uid} = upload;
    let newBlobHash = Promise.resolve(ZERO_HASH);
    let newBlobSize = 0;

    // we have already sent up to current_byte_position
    // we will send up to next_byte_position
    let next_byte_position = current_byte_position + next_blob_size;
    let blob_data = file.slice(current_byte_position, next_byte_position);

    // If there is still more file to send, estimate the new blob segment's details.
    if (current_byte_position < file.size) {
      let avgSpeed = !upload_speeds.length
        ? this.minBlobSize / XTR_TIME
        : upload_speeds.reduce((sum, t) => sum + t, 0) / upload_speeds.length;

      // we must also send the hash of the next blob,
      // so we must compute its size, between next_byte_position (the end
      // of this blob) and final_byte_position (the end of the next blob)
      newBlobSize = Math.min(
        Math.max(
          // we want at least the MIN_BLOB_SIZE
          this.minBlobSize,
          // but optimistically, based on the average speed
          Math.floor(XTR_TIME * avgSpeed)
        ),
        // but we'll stop at most at here
        this.maxBlobSize,
        // in fact we should stop when we hit the end of the file
        file.size - next_byte_position
      );

      let final_byte_position = next_byte_position + newBlobSize;
      let newBlob = file.slice(next_byte_position, final_byte_position);
      newBlobHash = hashBlob(newBlob);
    }

    return newBlobHash.then(newBlobHash => {
      let blobRequest = {

        action: 'blob',
        blob_data,
        next_blob_size: newBlobSize,
        next_blob_hash: newBlobHash,
        current_byte_position,
      }

      if (metis_uid) blobRequest[metis_uid] = metis_uid;
      
      return blobRequest;
    });
  }

  updateUpload(upload, synchronize = true) {
    const key = fileKey(upload);
    this.uploads = {...this.uploads, [key]: {...(this.uploads[key] || {}), ...upload}};

    if (synchronize) {
      this.synchronize();
    }

    return this.uploads[key];
  }

  removeUpload(upload) {
    const key = fileKey(upload);
    this.uploads = {...this.uploads};
    delete this.uploads[key];
    this.synchronize();
  }

  addUpload({upload}) {
    const key = fileKey(upload);
    this.numUploadFailuresByKey[key] = 0;
    return this.updateUpload({...upload, status: 'queued'});
  }

  cancelUpload({upload}) {
    const {url, project_name, file_name} = upload;
    this.removeUpload(upload);

    if (upload.status !== 'complete') {
      postUploadCancel(url, {project_name, file_name})
        .catch(err => this.reportError('Upload cancel failed', err));
    }
  }

  pauseUpload({upload}) {
    return this.updateUpload({...upload, status: 'paused'});
  }
}
