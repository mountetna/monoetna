import { setupWorker } from '../worker';

describe('setupWorker', () => {
  it('defines a dispatch method to postMessages', () => {
    const mockExecute = jest.fn();
    const mockWorker = {
      postMessage: jest.fn(),
      addEventListener: jest.fn()
    };

    expect(mockWorker.dispatch).toEqual(undefined);
    const worker = setupWorker(mockWorker, mockExecute);

    expect(worker).toEqual(mockWorker);
    expect(mockWorker.dispatch).not.toEqual(undefined);
    worker.dispatch({ foo: 'bar' });
    expect(mockWorker.postMessage).toHaveBeenCalledWith({ foo: 'bar' });
  });

  it('defines an error method to dispatch error messages', () => {
    const mockExecute = jest.fn();
    const mockWorker = {
      postMessage: jest.fn(),
      addEventListener: jest.fn()
    };

    expect(mockWorker.error).toEqual(undefined);
    const worker = setupWorker(mockWorker, mockExecute);

    expect(worker).toEqual(mockWorker);
    expect(mockWorker.error).not.toEqual(undefined);
    worker.error('It broke!');
    expect(mockWorker.postMessage).toHaveBeenCalledWith({
      type: 'WORKER_ERROR',
      worker: 'upload',
      message: 'It broke!'
    });
  });

  // It would be nice if we could attach a real event listenger, but
  //   then I think we would have to pass in a component for testing,
  //   not just a hash.
  it('attaches an event listener for messages', () => {
    const mockExecute = jest.fn();
    const mockWorker = {
      postMessage: jest.fn(),
      addEventListener: jest.fn()
    };

    const worker = setupWorker(mockWorker, mockExecute);

    expect(worker).toEqual(mockWorker);

    expect(mockWorker.addEventListener).toHaveBeenCalled();
  });
});
