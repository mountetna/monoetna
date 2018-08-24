import { createWorker } from '../workers';

let WORKERS = {
  upload: () => createWorker( require.resolve('../workers/uploader'))
}

// a middleware that dispatches commands to worker threads
const workRouter = () => {
  let workers = {
    upload: {}
  };
  return store => next => action => {
    let { type, work_type, name, ...args } = action;

    if (type != 'WORK') return next(action);

    if (!work_type in WORKERS) return;

    if (!(name in workers[work_type])) {
      let worker = WORKERS[work_type]();
      worker.addEventListener('message', ({data}) => store.dispatch(data))
      workers[work_type][name] = worker;
    }

    workers[work_type][name].postMessage(args);
  }
}

export default workRouter;
