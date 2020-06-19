import { stubUrl, cleanStubs } from '../../../spec/helpers';
import { SHOW_DIALOG } from '../../actions/message_actions';

import uploader from '../uploader';

describe('Uploader', () => {
  afterEach(() => {
    cleanStubs();
  });

  it('unqueue command puts some active uploads into the queue', () => {
    const mockWorker = {
      postMessage: jest.fn(),
      addEventListener: jest.fn()
    };

    uploader(mockWorker);

    const uploads = [
      { name: '1', status: 'active' },
      { name: '2', status: 'active' },
      { name: '3', status: 'active' },
      { name: '4', status: 'active' }
    ];
    mockWorker.unqueue({ uploads });

    expect(mockWorker.postMessage).toHaveBeenCalledWith({
      type: 'UPLOAD_STATUS',
      upload: { name: '4', status: 'active' },
      status: 'queued'
    });
  });

  it('unqueue command removes queued uploads and starts / continues them', () => {
    const mockWorker = {
      postMessage: jest.fn(),
      addEventListener: jest.fn()
    };

    uploader(mockWorker);

    // We set ``continue()`` to a mock here, because we'll
    //    test that method's functionality in separate specs.
    mockWorker.continue = jest.fn();

    const uploads = [
      { name: '1', status: 'active' },
      { name: '2', status: 'active' },
      { name: '3', status: 'queued' },
      { name: '4', status: 'queued' }
    ];
    mockWorker.unqueue({ uploads });

    expect(mockWorker.postMessage).toHaveBeenCalledWith({
      type: 'UPLOAD_STATUS',
      upload: { name: '3', status: 'queued' },
      status: 'active'
    });
    expect(mockWorker.continue).toHaveBeenCalledWith({
      upload: { name: '3', status: 'queued' }
    });
  });

  it('unqueue command does nothing when MAX_UPLOADS already active even with queued uploads', () => {
    const mockWorker = {
      postMessage: jest.fn(),
      addEventListener: jest.fn()
    };

    uploader(mockWorker);

    const uploads = [
      { name: '1', status: 'active' },
      { name: '2', status: 'active' },
      { name: '3', status: 'active' },
      { name: '4', status: 'queued' }
    ];
    mockWorker.unqueue({ uploads });

    expect(mockWorker.postMessage).not.toHaveBeenCalled();
  });

  it('unqueue command does nothing when no queued uploads and active <= MAX_UPLOADS', () => {
    const mockWorker = {
      postMessage: jest.fn(),
      addEventListener: jest.fn()
    };

    uploader(mockWorker);

    const uploads = [
      { name: '1', status: 'active' },
      { name: '2', status: 'active' },
      { name: '3', status: 'active' }
    ];
    mockWorker.unqueue({ uploads });

    expect(mockWorker.postMessage).not.toHaveBeenCalled();
  });

  it('start command sends POST request and dispatches message on success', (done) => {
    const mockWorker = {
      postMessage: jest.fn(),
      addEventListener: jest.fn()
    };

    uploader(mockWorker);

    const file = new File(['Once upon a time...'], 'legends.txt');

    const upload = {
      status: 'active',
      url: 'http://localhost',
      project_name: 'labors',
      file,
      file_size: file.size,
      file_name: file.name
    };

    global.fetch = require('node-fetch');

    stubUrl({
      verb: 'post',
      path: '/',
      request: {
        file_size: 19,
        next_blob_size: 19,
        next_blob_hash: '4aba72c5db1e9a47b575fdbe12a273de',
        action: 'start'
      },
      response: { upload },
      host: 'http://localhost'
    });

    mockWorker.start({ upload }).then(() => {
      expect(mockWorker.postMessage).toHaveBeenNthCalledWith(1, {
        type: 'UPLOAD_STATUS',
        upload: {
          upload: {
            status: 'active',
            url: 'http://localhost',
            project_name: 'labors',
            file: {}, // Can't expect a File object back from the server
            file_size: file.size,
            file_name: file.name
          }
        },
        status: 'queued'
      });

      expect(mockWorker.postMessage).toHaveBeenNthCalledWith(2, {
        type: 'UNQUEUE_UPLOADS'
      });
      done();
    });
  });

  it('continue command sends 0 blob request if current byte > file.size', (done) => {
    const mockWorker = {
      postMessage: jest.fn(),
      addEventListener: jest.fn()
    };

    uploader(mockWorker);

    const file = new File(['Once upon a time...'], 'legends.txt');

    const upload = {
      status: 'active',
      url: 'http://localhost',
      project_name: 'labors',
      file,
      file_size: file.size,
      file_name: file.name,
      current_byte_position: 50,
      next_blob_size: 19,
      upload_speeds: []
    };

    const responseUpload = {
      status: 'active',
      url: 'http://localhost',
      project_name: 'labors',
      file: {},
      current_byte_position: 50,
      next_blob_size: 0,
      upload_speeds: []
    };

    global.fetch = require('node-fetch');

    stubUrl({
      verb: 'post',
      path: '/',
      request: (body) => {
        // No great way to match on form data
        return body === '[object FormData]';
      },
      response: responseUpload,
      host: 'http://localhost'
    });

    mockWorker.continue({ upload }).then(() => {
      expect(mockWorker.postMessage).toHaveBeenNthCalledWith(1, {
        type: 'UPLOAD_SPEED',
        upload,
        upload_speed: 0
      });

      expect(mockWorker.postMessage).toHaveBeenNthCalledWith(2, {
        status: 'complete',
        type: 'UPLOAD_STATUS',
        upload: {
          ...upload,
          ...responseUpload
        }
      });

      expect(mockWorker.postMessage).toHaveBeenNthCalledWith(3, {
        type: 'UPLOAD_FILE_COMPLETED',
        upload: {
          ...upload,
          ...responseUpload
        }
      });
      done();
    });
  });

  it('continue command sends blob request if current byte < file.size', (done) => {
    const mockWorker = {
      postMessage: jest.fn(),
      addEventListener: jest.fn()
    };

    uploader(mockWorker);

    const file = new File(['Once upon a time...'], 'legends.txt');

    const upload = {
      status: 'active',
      url: 'http://localhost',
      project_name: 'labors',
      file,
      file_size: file.size,
      file_name: file.name,
      current_byte_position: 0,
      next_blob_size: 19,
      upload_speeds: []
    };

    const responseUpload = {
      status: 'active',
      url: 'http://localhost',
      project_name: 'labors',
      file: {},
      current_byte_position: 50,
      next_blob_size: 0,
      upload_speeds: []
    };

    global.fetch = require('node-fetch');

    stubUrl({
      verb: 'post',
      path: '/',
      request: (body) => {
        // No great way to match on form data
        return body === '[object FormData]';
      },
      response: responseUpload,
      host: 'http://localhost'
    });

    mockWorker.continue({ upload }).then(() => {
      // Can't guarantee the upload speed for this message, so
      //   we just check that it is called 3 times? Can we make
      //   this better?
      expect(mockWorker.postMessage).toHaveBeenCalledTimes(3);

      done();
    });
  });

  it('start command sends POST request and dispatches message on error', (done) => {
    const mockWorker = {
      postMessage: jest.fn(),
      addEventListener: jest.fn()
    };

    uploader(mockWorker);

    const file = new File(['Once upon a time...'], 'legends.txt');

    const upload = {
      status: 'active',
      url: 'http://localhost',
      project_name: 'labors',
      file,
      file_size: file.size,
      file_name: file.name
    };

    global.fetch = require('node-fetch');

    stubUrl({
      verb: 'post',
      path: '/',
      request: {
        file_size: 19,
        next_blob_size: 19,
        next_blob_hash: '4aba72c5db1e9a47b575fdbe12a273de',
        action: 'start'
      },
      status: 422,
      response: { error: 'Bad request' },
      host: 'http://localhost'
    });

    mockWorker.start({ upload }).then(() => {
      expect(mockWorker.postMessage).toHaveBeenNthCalledWith(1, {
        type: SHOW_DIALOG,
        dialog: {
          type: 'message',
          title: 'Upload failed',
          message: 'Bad request',
          message_type: 'warning'
        }
      });
      done();
    });
  });

  it('cancel command sends POST request and dispatches message on success', (done) => {
    const mockWorker = {
      postMessage: jest.fn(),
      addEventListener: jest.fn()
    };

    uploader(mockWorker);

    const upload = {
      status: 'active',
      url: 'http://localhost',
      project_name: 'labors',
      file_name: 'hydra.txt'
    };

    global.fetch = require('node-fetch');

    stubUrl({
      verb: 'post',
      path: '/',
      request: {
        project_name: 'labors',
        file_name: 'hydra.txt',
        action: 'cancel'
      },
      response: { success: true },
      host: 'http://localhost'
    });

    mockWorker.cancel({ upload }).then(() => {
      expect(mockWorker.postMessage).toHaveBeenCalledWith({
        type: 'UPLOAD_FILE_CANCELED',
        upload
      });
      done();
    });
  });

  it('cancel command sends POST request and dispatches message on error', (done) => {
    const mockWorker = {
      postMessage: jest.fn(),
      addEventListener: jest.fn()
    };

    uploader(mockWorker);

    const upload = {
      status: 'active',
      url: 'http://localhost',
      project_name: 'labors',
      file_name: 'hydra.txt'
    };

    global.fetch = require('node-fetch');

    stubUrl({
      verb: 'post',
      path: '/',
      request: {
        project_name: 'labors',
        file_name: 'hydra.txt',
        action: 'cancel'
      },
      status: 422,
      response: { error: true },
      host: 'http://localhost'
    });

    mockWorker.cancel({ upload }).then(() => {
      expect(mockWorker.postMessage).toHaveBeenCalledWith({
        dialog: {
          message: true,
          message_type: 'error',
          title: 'Upload cancel failed',
          type: 'message'
        },
        type: 'SHOW_DIALOG'
      });
      done();
    });
  });

  // Not sure how to test the setupWorker() execute callback,
  //   unless we export it somehow.
  xit('logs an exception if non-approved command is sent', () => {
    // Should this really throw an exception instead of just logging??
    global.console = jest.fn();
  });
});
