import React from 'react';
import { mockStore } from '../../spec/helpers';
import workDispatcher from '../work-dispatcher';
import { WORKERS } from '../work-dispatcher';
import { WORK, WORK_FAILED } from '../../upload/actions/upload_actions';

describe('workDispatcher', () => {
  let store;

  beforeEach(() => {
    store = mockStore({
      location: {
        path: '/labors/browse/monster/Nemean Lion'
      }
    });
  });

  it('ignores non-WORK actions', () => {
    const mockNext = jest.fn();

    const action = {
      type: 'NAP',
      work_type: 'Zzz'
    };
    workDispatcher()(store)(mockNext)(action);
    expect(mockNext).toHaveBeenCalledWith(action);
  });

  it('returns with non-upload WORK actions', () => {
    const mockNext = jest.fn();

    const action = {
      type: WORK,
      work_type: 'skip'
    };
    workDispatcher()(store)(mockNext)(action);
    expect(mockNext).not.toHaveBeenCalledWith(action);
  });

  it('posts message and creates new upload worker', () => {
    const mockNext = jest.fn();

    const otherArgs = {
      random: 'text'
    };

    const realCreateWorker = WORKERS['upload'];

    const mockWorker = {
      postMessage: jest.fn()
    };
    WORKERS['upload'] = jest.fn((dispatch) => {
      return mockWorker;
    });

    const action = {
      type: WORK,
      work_type: 'upload',
      command: otherArgs
    };
    workDispatcher()(store)(mockNext)(action);

    expect(mockWorker.postMessage).toHaveBeenCalledWith(otherArgs);

    WORKERS['upload'] = realCreateWorker;
  });
});
