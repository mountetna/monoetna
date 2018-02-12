/*
 * This class will initialize an upload to the Metis server from a web client.
 */

import md5 from 'md5';
import { startUpload } from '../actions/upload_actions';

module.exports = function(self) {
  let timeouts = 0; // The number of times an upload has timed out.
  let maxTimeouts = 5; // The number of attepts to upload a blob.
  const INITIAL_BLOB_SIZE = Math.pow(2, 10); // in bytes

  let dispatch = (action) => self.postMessage(action);

  // report errors
  let error = (message) => dispatch(
    { type: 'WORKER_ERROR', message }
  );

  /*
   * For message passing into the worker.
   */
  self.addEventListener('message', (event) => {
    let { file, url } = event.data;

    // Hash the next blob for security and error checking.
    let blob = file.slice(0, INITIAL_BLOB_SIZE);
    let fileReader = new FileReader();
    fileReader.onload = (event) => {
      let next_blob_hash = md5(fileReader.result);

      let upload = {
        next_blob_size: INITIAL_BLOB_SIZE,
        next_blob_hash
      };

      dispatch({ type: 'START_UPLOAD', upload, url});
    }
    fileReader.readAsArrayBuffer(blob);
  });
}
