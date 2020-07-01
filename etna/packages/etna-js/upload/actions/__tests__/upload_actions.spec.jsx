import 'core-js/stable';

import * as actions from '../upload_actions';
import { SHOW_DIALOG } from '../message_actions';

import { mockStore, stubUrl } from '../../../spec/helpers';

describe('upload actions', () => {
  it('unqueueUploads generates actions to unqueue all uploads from web worker', () => {
    const store = mockStore({});

    expect(store.getActions()).toEqual([]);

    actions.unqueueUploads()(store.dispatch, () => {
      return {
        directory: {
          uploads: [
            {
              file: 'upload-file.txt'
            }
          ]
        }
      };
    });

    expect(store.getActions()).toEqual([
      {
        type: actions.WORK,
        work_type: 'upload',
        command: 'unqueue',
        uploads: [
          {
            file: 'upload-file.txt'
          }
        ]
      }
    ]);
  });

  it('uploadFileCanceled generates action to cancel specific file upload', () => {
    const store = mockStore({});

    expect(store.getActions()).toEqual([]);

    const upload = {
      file: 'upload-file.txt'
    };

    actions.uploadFileCanceled({ upload })(store.dispatch);

    expect(store.getActions()).toEqual([
      {
        type: actions.REMOVE_UPLOAD,
        upload
      },
      {
        type: actions.UNQUEUE_UPLOADS
      }
    ]);
  });

  it('cancelUpload generates action to remove a completed file upload when cancelUpload called', () => {
    const store = mockStore({});

    expect(store.getActions()).toEqual([]);

    const upload = {
      file: 'upload-file.txt',
      status: 'complete'
    };

    actions.cancelUpload({ upload })(store.dispatch);

    expect(store.getActions()).toEqual([
      {
        type: actions.REMOVE_UPLOAD,
        upload
      }
    ]);
  });

  it('cancelUpload does nothing if user cancels the confirmation window', () => {
    window.confirm = jest.fn(() => false);
    const store = mockStore({});

    expect(store.getActions()).toEqual([]);

    const upload = {
      file: 'upload-file.txt',
      status: 'active'
    };

    actions.cancelUpload({ upload })(store.dispatch);

    expect(store.getActions()).toEqual([]);
  });

  it('cancelUpload generates work action if user confirms the cancellation', () => {
    window.confirm = jest.fn(() => true);
    const store = mockStore({});

    expect(store.getActions()).toEqual([]);

    const upload = {
      file: 'upload-file.txt',
      status: 'active'
    };

    actions.cancelUpload({ upload })(store.dispatch);

    expect(store.getActions()).toEqual([
      { type: actions.WORK, work_type: 'upload', command: 'cancel', upload }
    ]);
  });

  it('pauseUpload generates paused action if user pauses the upload', () => {
    const store = mockStore({});

    expect(store.getActions()).toEqual([]);

    const upload = {
      file: 'upload-file.txt',
      status: 'active'
    };

    actions.pauseUpload({ upload })(store.dispatch);

    expect(store.getActions()).toEqual([
      { type: actions.UPLOAD_STATUS, status: 'paused', upload }
    ]);
  });

  it('continueUpload generates actions to continue an upload', () => {
    const store = mockStore({});

    expect(store.getActions()).toEqual([]);

    const upload = {
      file: 'upload-file.txt',
      status: 'active'
    };

    actions.continueUpload({ upload })(store.dispatch);

    expect(store.getActions()).toEqual([
      { type: actions.UPLOAD_STATUS, status: 'active', upload },
      {
        type: actions.WORK,
        work_type: 'upload',
        command: 'continue',
        upload
      }
    ]);
  });

  it('uploadFileCompleted generates actions when file upload completed', () => {
    const store = mockStore({});

    expect(store.getActions()).toEqual([]);

    const upload = {
      file: 'upload-file.txt',
      status: 'active'
    };

    actions.uploadFileCompleted({ upload })(store.dispatch);

    expect(store.getActions()).toEqual([
      { type: actions.ADD_FILES, files: [upload.file] },
      {
        type: actions.UNQUEUE_UPLOADS
      }
    ]);
  });

  it('uploadBlobCompleted does nothing with a blob if the upload is not active', () => {
    const store = mockStore({});

    expect(store.getActions()).toEqual([]);

    global.CONFIG = {
      project_name: 'labors'
    };
    const file_name = 'hydra.txt';

    const getState = () => {
      return {
        directory: {
          uploads: {
            'labors:hydra.txt': {
              file: 'hydra.txt',
              status: 'paused'
            }
          }
        }
      };
    };

    actions.uploadBlobCompleted({ file_name })(store.dispatch, getState);

    expect(store.getActions()).toEqual([]);
  });

  it('uploadBlobCompleted dispatches work action to upload blob if the upload is active', () => {
    const store = mockStore({});

    expect(store.getActions()).toEqual([]);

    global.CONFIG = {
      project_name: 'labors'
    };
    const file_name = 'hydra.txt';

    const upload = {
      file: file_name,
      status: 'active'
    };

    const getState = () => {
      return {
        directory: {
          uploads: {
            'labors:hydra.txt': {
              file: file_name,
              status: 'active'
            }
          }
        }
      };
    };

    actions.uploadBlobCompleted({ file_name })(store.dispatch, getState);

    expect(store.getActions()).toEqual([
      {
        type: actions.WORK,
        work_type: 'upload',
        command: 'continue',
        upload
      }
    ]);
  });

  it('uploadStarted does not start an upload if the upload is not active', () => {
    const store = mockStore({});

    expect(store.getActions()).toEqual([]);

    global.CONFIG = {
      project_name: 'labors'
    };
    const file_name = 'hydra.txt';

    const getState = () => {
      return {
        directory: {
          uploads: {
            'labors:hydra.txt': {
              file: 'hydra.txt',
              status: 'paused'
            }
          }
        }
      };
    };

    actions.uploadStarted({ file_name })(store.dispatch, getState);

    expect(store.getActions()).toEqual([]);
  });

  it('uploadStarted dispatches work action to start upload if the upload is active', () => {
    const store = mockStore({});

    expect(store.getActions()).toEqual([]);

    global.CONFIG = {
      project_name: 'labors'
    };
    const file_name = 'hydra.txt';

    const upload = {
      file: file_name,
      status: 'active'
    };

    const getState = () => {
      return {
        directory: {
          uploads: {
            'labors:hydra.txt': upload
          }
        }
      };
    };

    actions.uploadStarted({ file_name })(store.dispatch, getState);

    expect(store.getActions()).toEqual([
      {
        type: actions.WORK,
        work_type: 'upload',
        command: 'continue',
        upload
      }
    ]);
  });

  it('fileSelected dispatches work actions to start upload', async () => {
    const store = mockStore({});

    expect(store.getActions()).toEqual([]);

    const file_name = 'hydra.txt';

    const upload = {
      file: file_name,
      status: 'active'
    };

    const getState = () => {
      return {
        directory: {
          uploads: {
            'labors:profiles/hydra.txt': upload
          }
        }
      };
    };

    const url = 'https://www.mountolympus.org';

    stubUrl({
      verb: 'post',
      path: '/authorize/upload',
      request: {
        project_name: 'labors',
        file_path: 'profiles/hydra.txt',
        bucket_name: 'monsters'
      },
      response: { url },
      host: 'http://localhost'
    });

    const dispatch = jest.fn();

    await actions.fileSelected({
      file: {
        name: file_name
      },
      folder_name: 'profiles',
      bucket_name: 'monsters'
    })(dispatch, getState);

    expect(dispatch).toHaveBeenNthCalledWith(1, {
      type: actions.ADD_UPLOAD,
      project_name: CONFIG.project_name,
      file: {
        name: file_name
      },
      file_name: 'profiles/hydra.txt',
      url
    });

    expect(dispatch).toHaveBeenNthCalledWith(2, {
      type: actions.WORK,
      work_type: 'upload',
      command: 'start',
      upload
    });
  });

  // I can't figure out how to test the second sad path, because
  //   of the promise handling in message_actions.errorMessage...
  it('fileSelected throws error message if authorization fails', (done) => {
    const store = mockStore({});

    expect(store.getActions()).toEqual([]);

    const file_name = 'hydra.txt';

    const upload = {
      file: file_name,
      status: 'active'
    };

    const getState = () => {
      return {
        directory: {
          uploads: {
            'labors:profiles/hydra.txt': upload
          }
        }
      };
    };

    stubUrl({
      verb: 'post',
      path: '/authorize/upload',
      request: {
        project_name: 'labors',
        file_path: 'profiles/hydra.txt',
        bucket_name: 'monsters'
      },
      status: 422,
      response: { error: 'Error!' },
      host: 'http://localhost'
    });

    const dispatch = jest.fn();

    actions
      .fileSelected({
        file: {
          name: file_name
        },
        folder_name: 'profiles',
        bucket_name: 'monsters'
      })(dispatch, getState)
      .then(() => {
        expect(dispatch).toHaveBeenNthCalledWith(1, {
          type: SHOW_DIALOG,
          dialog: {
            type: 'message',
            title: 'Upload failed',
            message: 'Error!',
            message_type: 'warning'
          }
        });
        done();
      });
  });
});
