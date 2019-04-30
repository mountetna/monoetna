export const setupWorker = function(self, commands) {
  let worker = {
    dispatch: (action) => self.postMessage(action),
    error: (message) => worker.dispatch(
      { type: 'WORKER_ERROR', worker: 'upload', message }
    )
  };

  self.addEventListener('message', ({data}) => {
    let { command } = data;

    // Check that the incoming data has a valid command.
    if (command in commands)
      commands[command](data);
    else {
      worker.error(`Invalid command ${command}`);
    }
  });
  return worker;
}
