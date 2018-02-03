import * as Redux from 'redux';
import * as ReduxLogger from 'redux-logger';
import asyncRouter from './async-router';

const DEFAULT_STATE = {
  fileData: {
    fileList: [],
    fileUploads: [],
    fileFails: []
  },

  userInfo: {
    userId: null,
    userEmail: '',
    authToken: '',
    firstName: '',
    lastName: '',
    permissions: [],

    masterPerms: false,

    loginStatus: false,
    loginError: false,
    loginErrorMsg: 'Invalid sign in.'
  }
};

const createStore = (reducers, actions) => {
  let reducer = Redux.combineReducers(reducers);

  let middleWares = [
    asyncRouter(actions)
  ];

  if(process.env.NODE_ENV != 'production') middleWares.push(ReduxLogger.createLogger());

  let store = Redux.applyMiddleware(...middleWares)(Redux.createStore)(reducer, DEFAULT_STATE);

  return store;
}

export default createStore;
