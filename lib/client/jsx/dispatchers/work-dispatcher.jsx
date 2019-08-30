import { createWorker } from '../workers';

let WORKERS = {
  upload: () => createWorker( require.resolve('../workers/uploader'))
}

// a middleware that dispatches commands to worker threads
const workDispatcher = () => {
  let workers = {
    upload: {}
  };
  return store => next => action => {
    let { type, work_type, worker_name, ...args } = action;

    if (type != 'WORK') return next(action);

    if (!work_type in WORKERS) return;

    if (!(worker_name in workers[work_type])) {
      let worker = WORKERS[work_type]();
      worker.addEventListener('message', ({data}) => store.dispatch(data))
      workers[work_type][worker_name] = worker;
    }

    workers[work_type][worker_name].postMessage(args);
  }
}

export default workDispatcher;
