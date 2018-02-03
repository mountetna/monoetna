export const queueUpload = () => (dispatch) => {
  this.startUpload();
}

export const pauseUpload = () => (dispatch) => {
  this.uploadWorker.postMessage({ command: 'pause' });
}

export const cancelUpload = () => (dispatch) => {
  if(!confirm('Are you sure you want to remove this upload?')) return;
  this.uploadWorker.postMessage({ command: 'cancel' });
}

export const recoverUpload = () => (dispatch) => {
  this.recoverUpload(action.uploadFile, action.fileMetadata);
}
