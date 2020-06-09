export const setupWorker = function(worker, execute) {
    Object.assign(worker, {
      dispatch: (action) => worker.postMessage(action),
      error: (message) => worker.dispatch(
        { type: 'WORKER_ERROR', worker: 'upload', message }
      )
    });

    worker.addEventListener('message', ({data}) => execute(data));
    return worker;
  }
