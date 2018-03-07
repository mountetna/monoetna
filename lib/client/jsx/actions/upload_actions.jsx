
export const queueUpload = ({upload}) => (dispatch) => {
  dispatch({ type: 'WORK', worker: 'upload', command: 'start', upload });
}

export const pauseUpload = ({upload}) => (dispatch) => {
  dispatch({ type: 'WORK', worker: 'upload', command: 'pause', upload });
}

export const cancelUpload = ({upload}) => (dispatch) => {
  if(!confirm('Are you sure you want to remove this upload?')) return;
  dispatch({ type: 'WORK', worker: 'upload', command: 'cancel', upload });
}
