import * as Redux from 'redux';
import * as ReduxLogger from 'redux-logger';
import fileData from './metis-reducer';
import JanusLogReducer from './janus-log-reducer';
import LastActionReducer from './last-action-reducer';

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

const createStore = () => {
  let janusLogReducer = new JanusLogReducer();
  let lastAction = new LastActionReducer();

  let reducer = Redux.combineReducers({
    fileData,
    userInfo: janusLogReducer.reducer(),
    lastAction: lastAction.reducer()
  });


  let middleWares = []
  if(process.env.NODE_ENV != 'production') middleWares.push(ReduxLogger.createLogger());

  let store = Redux.applyMiddleware(...middleWares)(Redux.createStore)(reducer, DEFAULT_STATE);

  return store;
}

export default createStore;
