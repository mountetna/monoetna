import {createWorker, terminateWorker} from '../web-workers/index';
import { WORK, WORK_FAILED } from '../upload/actions/upload_actions';
import {appSubscription} from "../utils/subscription";

// Export for testing
export const WORKERS = {
  upload: (dispatch) => createWorker( dispatch, require.resolve('../upload/workers/uploader'))
}

// a middleware that dispatches commands to worker threads
const workDispatcher = () => {
  let active_workers = {};

  appSubscription.addCleanup(function clearAllWorkers() {
    Object.values(active_workers).forEach(terminateWorker)
  });

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

export default workDispatcher;
