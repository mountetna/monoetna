import { createWorker } from '../workers/index';
import { WORK, WORK_FAILED } from '../actions/upload_actions';

let WORKERS = {
  upload: (dispatch) => createWorker( dispatch, require.resolve('../workers/uploader'))
}

// a middleware that dispatches commands to worker threads
const workDispatcher = () => {
  let active_workers = { };

  return store => next => action => {
    let { type, work_type, ...args } = action;

    if (type != WORK) return next(action);

    if (!work_type in WORKERS) return;

    if (!(work_type in active_workers)) {
      try {
        let worker = WORKERS[work_type](store.dispatch);
        active_workers[work_type] = worker;
      } catch(e) {
        store.dispatch({
          type: WORK_FAILED, work_type, ...args
        });
        return;
      }
    }

    active_workers[work_type].postMessage(args);
  }
}

export default workDispatcher;
