// a middleware that dispatches commands to worker threads

const workRouter = workers => {
  return store => {
    Object.values(workers).forEach(
      worker=>worker.addEventListener('message', ({data}) => store.dispatch(data))
    );
    return next => action => {
      let { type, worker, ...args } = action;

      if (type != 'WORK') return next(action);

      if (workers[worker]) workers[worker].postMessage(args);
    }
  }
}

export default workRouter;
