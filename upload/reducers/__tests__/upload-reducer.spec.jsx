import React from 'react';
import { mockStore } from '../../../spec/helpers';
import renderer from 'react-test-renderer';
import UploadReducer from '../upload-reducer';

import {
  UPLOAD_SPEED,
  UPLOAD_STATUS,
  ADD_UPLOAD,
  REMOVE_UPLOAD
} from '../../actions/upload_actions';

describe('Upload Reducer', () => {
  let store;

  beforeEach(() => {
    store = mockStore({
      location: {
        path: '/labors/browse/monster/Nemean Lion'
      }
    });
  });

  it('returns old uploads for unrecognized action', () => {
    const oldUploads = {
      foo: 'bar'
    };
    const result = UploadReducer(oldUploads, { type: 'unrecognized' });
    expect(result).toEqual(oldUploads);
  });

  it('updates the upload speeds for UPLOAD_SPEED action', () => {
    const oldUploads = {
      'labors:hydra.txt': {
        project_name: 'labors',
        file_name: 'hydra.txt'
      },
      'labors:lion.txt': {}
    };
    const action = {
      type: UPLOAD_SPEED,
      upload_speed: 500,
      upload: {
        project_name: 'labors',
        file_name: 'hydra.txt',
        upload_speeds: [100, 200, 300, 400]
      }
    };
    const result = UploadReducer(oldUploads, action);
    expect(result).toEqual({
      'labors:hydra.txt': {
        project_name: 'labors',
        file_name: 'hydra.txt',
        upload_speeds: [100, 200, 300, 400, 500]
      },
      'labors:lion.txt': {}
    });
  });

  it('updates the upload status for UPLOAD_STATUS action', () => {
    const oldUploads = {
      'labors:hydra.txt': {
        project_name: 'labors',
        file_name: 'hydra.txt',
        upload_speeds: [100, 200, 300, 400]
      },
      'labors:lion.txt': {}
    };
    const action = {
      type: UPLOAD_STATUS,
      status: 'active',
      upload: {
        project_name: 'labors',
        file_name: 'hydra.txt',
        current_byte_position: 1024,
        next_blob_size: 8096,
        next_blob_hash: '0x123'
      }
    };
    const result = UploadReducer(oldUploads, action);
    expect(result).toEqual({
      'labors:hydra.txt': {
        project_name: 'labors',
        file_name: 'hydra.txt',
        status: 'active',
        current_byte_position: 1024,
        next_blob_size: 8096,
        next_blob_hash: '0x123',
        upload_speeds: [100, 200, 300, 400]
      },
      'labors:lion.txt': {}
    });
  });

  it('adds a new upload for ADD_UPLOAD action', () => {
    const oldUploads = {
      'labors:hydra.txt': {
        project_name: 'labors',
        file_name: 'hydra.txt',
        upload_speeds: [100, 200, 300, 400]
      },
      'labors:lion.txt': {}
    };
    const action = {
      type: ADD_UPLOAD,
      project_name: 'labors',
      file_name: 'stables.txt',
      file: {
        size: 1024
      },
      url: 'https://www.mountolympus.org'
    };
    const result = UploadReducer(oldUploads, action);
    expect(result).toEqual({
      'labors:stables.txt': {
        project_name: 'labors',
        file_name: 'stables.txt',
        file: {
          size: 1024
        },
        url: 'https://www.mountolympus.org',
        file_size: 1024,
        current_byte_position: 0,
        status: 'preparing',
        upload_speeds: []
      },
      'labors:lion.txt': {},
      'labors:hydra.txt': {
        project_name: 'labors',
        file_name: 'hydra.txt',
        upload_speeds: [100, 200, 300, 400]
      }
    });
  });

  it('removes an upload for REMOVE_UPLOAD action', () => {
    const oldUploads = {
      'labors:hydra.txt': {
        project_name: 'labors',
        file_name: 'hydra.txt',
        upload_speeds: [100, 200, 300, 400]
      },
      'labors:lion.txt': {}
    };
    const action = {
      type: REMOVE_UPLOAD,
      upload: {
        project_name: 'labors',
        file_name: 'hydra.txt'
      }
    };
    const result = UploadReducer(oldUploads, action);
    expect(result).toEqual({
      'labors:lion.txt': {}
    });
  });
});
