import {createWorker, terminateWorker} from '../upload/workers/index';
import { WORK, WORK_FAILED } from '../upload/actions/upload_actions';

// Export for testing
export const WORKERS = {
  upload: (dispatch) => createWorker( dispatch, require.resolve('../upload/workers/uploader'))
}

// For tests to disable all active workers and prevent test polution of async operations.
const clearAllWorkersQueue = [];

// a middleware that dispatches commands to worker threads
const workDispatcher = () => {
  let active_workers = { };


  function clearAllWorkers() {
    Object.values(active_workers).forEach(terminateWorker)
  }

  clearAllWorkersQueue.push(clearAllWorkers);

  return store => next => action => {
    let { type, work_type, command } = action;

    if (type != WORK) return next(action);

    if (!work_type in WORKERS) return;

    if (!(work_type in active_workers)) {
      try {
        let worker = WORKERS[work_type](store.dispatch);
        active_workers[work_type] = worker;
      } catch(e) {
        store.dispatch({
          type: WORK_FAILED, work_type, command
        });
        return;
      }
    }

    active_workers[work_type].postMessage(command);
  }
}

export function terminateAllWorkers() {
  clearAllWorkersQueue.forEach(f => f())
  clearAllWorkersQueue.length = 0;
}

export default workDispatcher;
