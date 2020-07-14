import * as api from '../upload_api';

import { mockStore, stubUrl, cleanStubs } from '../../../spec/helpers';

describe('upload api', () => {
  it('postAuthorizeUpload makes a request to the authorize upload route', () => {
    stubUrl({
      verb: 'post',
      path: '/authorize/upload',
      request: {
        project_name: global.CONFIG.project_name,
        file_path: 'profiles/hydra.txt',
        bucket_name: 'monsters'
      },
      response: { success: true },
      host: 'http://localhost'
    });

    const response = api.postAuthorizeUpload(
      'http://localhost',
      'labors',
      'monsters',
      'profiles/hydra.txt'
    );

    return response.then((res) => {
      expect(res).toEqual({ success: true });
    });
  });

  it('postUploadStart makes a request to the upload route with start action', () => {
    stubUrl({
      verb: 'post',
      path: '/upload',
      request: {
        action: 'start',
        file: {
          name: 'hydra.txt'
        }
      },
      response: { success: true },
      host: 'http://localhost'
    });

    const response = api.postUploadStart('http://localhost/upload', {
      file: {
        name: 'hydra.txt'
      }
    });

    return response.then((res) => {
      expect(res).toEqual({ success: true });
    });
  });

  it('postUploadCancel makes a request to the upload route with cancel action', () => {
    stubUrl({
      verb: 'post',
      path: '/upload',
      request: {
        action: 'cancel',
        file: {
          name: 'hydra.txt'
        }
      },
      response: { success: true },
      host: 'http://localhost'
    });

    const response = api.postUploadCancel('http://localhost/upload', {
      file: {
        name: 'hydra.txt'
      }
    });

    return response.then((res) => {
      expect(res).toEqual({ success: true });
    });
  });

  it('postUploadBlob makes a request to the upload route with data blob', () => {
    stubUrl({
      verb: 'post',
      path: '/upload',
      request: /action\=blob/,
      response: { success: true },
      host: 'http://localhost'
    });

    const response = api.postUploadBlob('http://localhost/upload', {
      file: 'hydra.txt'
    });

    return response.then((res) => {
      expect(res).toEqual({ success: true });
    });
  });
});
