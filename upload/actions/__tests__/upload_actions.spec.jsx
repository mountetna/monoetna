import 'core-js/stable';

import * as actions from '../upload_actions';
import { SHOW_DIALOG } from '../message_actions';

import {delay, mockStore, stubUrl} from '../../../spec/helpers';
import {ADD_FILES, unpauseUpload, WORK} from "../upload_actions";
import {AddUploadCommand, Upload, UPLOAD_COMPLETE} from "../../workers/uploader";
import workDispatcher from "../../../dispatchers/work-dispatcher";
import asyncDispatcher from "../../../dispatchers/async-dispatcher";

describe('upload actions', () => {
  const emptyFile = new File([], 'empty.txt');
  const emptyUpload = Upload({ file_name: 'empty.txt', project_name: CONFIG.project_name, url: 'http://localhost/empty', file: emptyFile });

  describe('worker integration', () => {
    it('shows dialogs for upload failures', async () => {
      const { getActions, dispatch } = mockStore({}, [asyncDispatcher(actions), workDispatcher()]);

      // Test assumption: empty file uploads still pass through start and blob once.
      // This just makes the stubbing very simple.
      // Once for the start
      stubUrl({
        verb: 'post',
        path: '/empty',
        response: { error: 'Oh no' },
        status: 401,
        request: /.*/,
      })

      dispatch(unpauseUpload({ upload: emptyUpload }));

      for (let i = 0; i < 10; ++i) {
        await delay(100);
        if (getActions().map(({ type }) => type).includes(SHOW_DIALOG)) break;
      }

      expect(getActions()).toContainEqual({
        type: SHOW_DIALOG,
        dialog: {
          message: 'Oh no',
          message_type: 'error',
          title: 'Upload of file empty.txt',
          type: 'message'
        }
      });
    })

    it('upload completion adds files', async () => {
      const { getActions, dispatch } = mockStore({}, [asyncDispatcher(actions), workDispatcher()]);

      const { file_name, project_name, url, current_byte_position } = emptyUpload;

      // Test assumption: empty file uploads still pass through start and blob once.
      // This just makes the stubbing very simple.
      // Once for the start
      stubUrl({
        verb: 'post',
        path: '/empty',
        response: { file_name, project_name, url, current_byte_position },
        request: /.*/,
      })
      // Once for the blob
      stubUrl({
        verb: 'post',
        path: '/empty',
        response: { file_name, project_name, url, current_byte_position },
        request: /.*/,
      })

      dispatch(unpauseUpload({ upload: emptyUpload }));

      for (let i = 0; i < 10; ++i) {
        await delay(100);
        if (getActions().map(({ type }) => type).includes(UPLOAD_COMPLETE)) break;
      }

      expect(getActions()).toContainEqual({
        type: ADD_FILES,
        files: [emptyFile],
      });
    })
  })

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
      type: WORK,
      work_type: 'upload',
      command: AddUploadCommand(Upload({
        project_name: CONFIG.project_name,
        file: {
          name: file_name
        },
        file_name: 'profiles/hydra.txt',
        url
      }))
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
