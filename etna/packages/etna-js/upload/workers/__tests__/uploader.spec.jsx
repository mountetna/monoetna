import {stubUrl, cleanStubs} from '../../../spec/helpers';
import {SHOW_DIALOG} from '../../actions/message_actions';

import uploader from '../uploader';

class TestUploader {
  constructor() {
    this.postMessage = jest.fn();
    this.addEventListener = jest.fn();
    uploader(this);
  }
}

describe('Uploader', () => {
  afterEach(() => {
    cleanStubs();
  });

  it('unqueue command puts some active uploads into the queue', () => {
    const uploader = new TestUploader();

    const uploads = [
      {name: '1', status: 'active'},
      {name: '2', status: 'active'},
      {name: '3', status: 'active'},
      {name: '4', status: 'active'}
    ];
    uploader.unqueue({uploads});

    expect(uploader.postMessage).toHaveBeenCalledWith({
      type: 'UPLOAD_STATUS',
      upload: {name: '4', status: 'active'},
      status: 'queued'
    });
  });

  it('unqueue command removes queued uploads and starts / continues them', () => {
    const uploader = new TestUploader();

    // We set ``continue()`` to a mock here, because we'll
    //    test that method's functionality in separate specs.
    uploader.continue = jest.fn();

    const uploads = [
      {name: '1', status: 'active'},
      {name: '2', status: 'active'},
      {name: '3', status: 'queued'},
      {name: '4', status: 'queued'}
    ];
    uploader.unqueue({uploads});

    expect(uploader.postMessage).toHaveBeenCalledWith({
      type: 'UPLOAD_STATUS',
      upload: {name: '3', status: 'queued'},
      status: 'active'
    });
    expect(uploader.continue).toHaveBeenCalledWith({
      upload: {name: '3', status: 'queued'}
    });
  });

  it('unqueue command does nothing when MAX_UPLOADS already active even with queued uploads', () => {
    const uploader = new TestUploader();

    const uploads = [
      {name: '1', status: 'active'},
      {name: '2', status: 'active'},
      {name: '3', status: 'active'},
      {name: '4', status: 'queued'}
    ];
    uploader.unqueue({uploads});

    expect(uploader.postMessage).not.toHaveBeenCalled();
  });

  it('unqueue command does nothing when no queued uploads and active <= MAX_UPLOADS', () => {
    const uploader = new TestUploader();

    const uploads = [
      {name: '1', status: 'active'},
      {name: '2', status: 'active'},
      {name: '3', status: 'active'}
    ];
    uploader.unqueue({uploads});

    expect(uploader.postMessage).not.toHaveBeenCalled();
  });

  it('start command sends POST request and dispatches message on success', (done) => {
    const uploader = new TestUploader();

    const file = new File(['Once upon a time...'], 'legends.txt');

    const upload = {
      status: 'active',
      url: 'http://localhost',
      project_name: 'labors',
      file,
      file_size: file.size,
      file_name: file.name
    };

    stubUrl({
      verb: 'post',
      path: '/',
      request: {
        file_size: 19,
        next_blob_size: 19,
        next_blob_hash: '4aba72c5db1e9a47b575fdbe12a273de',
        action: 'start',
        reset: false,
      },
      response: {upload},
      host: 'http://localhost'
    });

    uploader.start({upload}).then(() => {
      expect(uploader.postMessage).toHaveBeenNthCalledWith(1, {
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

      expect(uploader.postMessage).toHaveBeenNthCalledWith(2, {
        type: 'UNQUEUE_UPLOADS'
      });
      done();
    });
  });

  it('continue command sends 0 blob request if current byte > file.size', (done) => {
    const uploader = new TestUploader();

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

    stubUrl({
      verb: 'post',
      path: '/',
      request: /action\=blob/,
      response: responseUpload,
      host: 'http://localhost'
    });

    uploader.continue({upload}).then(() => {
      expect(uploader.postMessage).toHaveBeenNthCalledWith(1, {
        type: 'UPLOAD_SPEED',
        upload,
        upload_speed: 0
      });

      expect(uploader.postMessage).toHaveBeenNthCalledWith(2, {
        status: 'complete',
        type: 'UPLOAD_STATUS',
        upload: {
          ...upload,
          ...responseUpload
        }
      });

      expect(uploader.postMessage).toHaveBeenNthCalledWith(3, {
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
    const uploader = new TestUploader();

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

    stubUrl({
      verb: 'post',
      path: '/',
      request: /action\=blob/,
      response: responseUpload,
      host: 'http://localhost'
    });

    uploader.continue({upload}).then(() => {
      // Can't guarantee the upload speed for this message, so
      //   we just check that it is called 3 times? Can we make
      //   this better?
      expect(uploader.postMessage).toHaveBeenCalledTimes(3);

      done();
    });
  });

  it('start command sends POST request and dispatches message on error', (done) => {
    const uploader = new TestUploader();

    const file = new File(['Once upon a time...'], 'legends.txt');

    const upload = {
      status: 'active',
      url: 'http://localhost',
      project_name: 'labors',
      file,
      file_size: file.size,
      file_name: file.name
    };

    stubUrl({
      verb: 'post',
      path: '/',
      request: {
        file_size: 19,
        next_blob_size: 19,
        next_blob_hash: '4aba72c5db1e9a47b575fdbe12a273de',
        action: 'start',
        reset: false,
      },
      status: 500,
      response: {error: 'Bad request'},
      host: 'http://localhost'
    });

    uploader.start({upload}).then(() => {
      expect(uploader.postMessage).toHaveBeenNthCalledWith(1, {
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
    const uploader = new TestUploader();

    const upload = {
      status: 'active',
      url: 'http://localhost',
      project_name: 'labors',
      file_name: 'hydra.txt'
    };

    stubUrl({
      verb: 'post',
      path: '/',
      request: {
        project_name: 'labors',
        file_name: 'hydra.txt',
        action: 'cancel'
      },
      response: {success: true},
      host: 'http://localhost'
    });

    uploader.cancel({upload}).then(() => {
      expect(uploader.postMessage).toHaveBeenCalledWith({
        type: 'UPLOAD_FILE_CANCELED',
        upload
      });
      done();
    });
  });

  it('cancel command sends POST request and dispatches message on error', (done) => {
    const uploader = new TestUploader();

    const upload = {
      status: 'active',
      url: 'http://localhost',
      project_name: 'labors',
      file_name: 'hydra.txt'
    };

    stubUrl({
      verb: 'post',
      path: '/',
      request: {
        project_name: 'labors',
        file_name: 'hydra.txt',
        action: 'cancel'
      },
      status: 422,
      response: {error: true},
      host: 'http://localhost'
    });

    uploader.cancel({upload}).then(() => {
      expect(uploader.postMessage).toHaveBeenCalledWith({
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

  describe('continue', () => {
    const file = new File(['A file upload contents'], 'myfile.txt');

    const upload = {
      status: 'active',
      url: 'http://localhost/upload-123',
      project_name: 'labors',
      file,
      file_size: file.size,
      file_name: file.name,
      current_byte_position: 0,
      next_blob_size: 10,
      upload_speeds: []
    };

    describe('on 422 response', () => {
      it('after maxTimeouts failures for the same url stalls', async () => {
        const uploader = new TestUploader();

        const mockUploadAttempt = () => stubUrl({
          verb: 'post',
          path: '/upload-123',
          status: 422,
          request: /action\=blob/,
          host: 'http://localhost',
          response: {},
        })

        const mockRestartAttempt = () => stubUrl({
          verb: 'post',
          path: '/upload-123',
          status: 200,
          request: (json) => json.reset && json.action === 'start',
          host: 'http://localhost',
          response: {upload},
        });

        for (let i = 0; i < uploader.maxTimeouts - 1; ++i) {
          mockUploadAttempt().then(mockRestartAttempt);
          await uploader.continue({ upload });

          expect(uploader.postMessage.mock.calls.map(([data]) => data.type)).toEqual([
            "UPLOAD_STATUS",
            "UNQUEUE_UPLOADS",
          ]);
          uploader.postMessage.mock.calls.length = 0;
        }

        mockUploadAttempt();
        await uploader.continue({ upload });
        expect(uploader.postMessage.mock.calls.map(([data]) => data.type)).toEqual([
          "UPLOAD_TIMEOUT",
        ]);
      });

      it('attempts a reset', async () => {
        const uploader = new TestUploader();
        const requestChain = stubUrl({
          verb: 'post',
          path: '/upload-123',
          status: 422,
          request: /action\=blob/,
          host: 'http://localhost',
          response: {},
        }).then(() => stubUrl({
          verb: 'post',
          path: '/upload-123',
          status: 200,
          request: (json) => json.reset && json.action === 'start',
          host: 'http://localhost',
          response: {upload},
        }));


        await uploader.continue({upload});
        await expect(requestChain).toResolve();

        // No blob completes; those would only occur after the unqueue retry is attempted.
        expect(uploader.postMessage.mock.calls.map(([data]) => data.type)).toEqual([
          "UPLOAD_STATUS",
          "UNQUEUE_UPLOADS",
        ]);
      });
    });
  });

  describe('timeout', () => {
    const upload1 = {url: 'b'};
    const upload2 = {url: 'c'};

    it('is false for any given upload url until called a max number of times', () => {
      const uploader = new TestUploader();

      expect(uploader.maxTimeouts).toBeGreaterThan(0);
      for (let i = 0; i < uploader.maxTimeouts - 1; ++i) {
        expect(uploader.timeout(upload1)).toBeFalsy();
      }

      expect(uploader.timeout(upload1)).toBeTruthy();
      expect(uploader.timeout(upload1)).toBeTruthy();
      expect(uploader.timeout(upload2)).toBeFalsy();

      uploader.reset();
      expect(uploader.timeout(upload1)).toBeFalsy();
    })
  });

  // Not sure how to test the setupWorker() execute callback,
  //   unless we export it somehow.
  xit('logs an exception if non-approved command is sent', () => {
    // Should this really throw an exception instead of just logging??
    global.console = jest.fn();
  });
});
