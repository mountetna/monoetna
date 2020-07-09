import work from 'webworkify-webpack'

const maxWorkers = window.navigator.hardwareConcurrency || 4;

// keep track of workers in use
let workerMap = new Map();

export function createWorker(dispatch, script) {
  if (workerMap.size >= maxWorkers) {
    console.warn(`too many workers max: ${maxWorkers} - make sure to clean them up with terminateWorker`);
  }

  let worker = work(script);
  worker.messageHandler = ({data}) => dispatch(data);
  worker.addEventListener('message', worker.messageHandler);
  workerMap.set(worker, true);
  return worker;
}

export function terminateWorker(worker) {
  console.log({ worker });
  worker.removeEventListener('message', worker.messageHandler);
  worker.terminate();
  workerMap.delete(worker);
}
