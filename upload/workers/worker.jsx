// Provides observable dispatch and subscribe functions for listening to messages and post replies.
export const setupWorker = function (worker) {
  const dispatch = (action) => worker.postMessage(action);
  const subscribe = (f) => worker.addEventListener('message', ({data}) => f(data));
  return {subscribe, dispatch};
}
