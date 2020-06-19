import * as actions from '../message_actions';

import { mockStore } from '../../spec/helpers';

describe('message actions', () => {
  it('generates action to show a dialog', () => {
    const action = actions.message(
      'ALERT',
      'This is a message dialog',
      'I hope you read this.'
    );

    expect(action).toEqual({
      type: actions.SHOW_DIALOG,
      dialog: {
        type: 'message',
        title: 'This is a message dialog',
        message: 'I hope you read this.',
        message_type: 'ALERT'
      }
    });
  });

  it('dispatches error message action to show a dialog with a Promise', () => {
    const store = mockStore({});

    expect(store.getActions()).toEqual([]);

    const messageHandler = (msg) => msg;

    const promiseResponse = new Promise((resolve, reject) => {
      resolve({ error: 'Another message.' });
    });

    const dispatch = (payload) => {
      expect(payload).toEqual({
        type: actions.SHOW_DIALOG,
        dialog: {
          type: 'message',
          title: 'An error message',
          message: 'Another message.',
          message_type: 'ERROR'
        }
      });
    };

    actions.errorMessage(
      dispatch,
      'ERROR',
      'An error message',
      messageHandler
    )(promiseResponse);
  });

  it('dispatches error message action to show a dialog without a Promise', () => {
    const store = mockStore({});

    expect(store.getActions()).toEqual([]);

    const messageHandler = (msg) => msg;

    actions.errorMessage(
      store.dispatch,
      'ERROR',
      'An error message',
      messageHandler
    )('Another message.');

    expect(store.getActions()).toEqual([
      {
        type: actions.SHOW_DIALOG,
        dialog: {
          type: 'message',
          title: 'An error message',
          message: 'Another message.',
          message_type: 'ERROR'
        }
      }
    ]);
  });
});
