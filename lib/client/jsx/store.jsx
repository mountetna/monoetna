import * as Redux from 'redux';
import * as ReduxLogger from 'redux-logger';
import directory from './reducers/directory-reducer';
import user from './reducers/janus-log-reducer';

import * as fileActions from './actions/file_actions';
import * as uploadActions from './actions/upload_actions';
import * as userActions from './actions/user_actions';

import { createWorker } from './workers';
import asyncRouter from './routers/async-router';
import workRouter from './routers/work-router';

const createStore = () => {
  let reducers = {
    directory,
    user
  };

  // action handlers to import
  let actions = {
    ...fileActions,
    ...uploadActions,
    ...userActions

    // here you may define aliases to other actions,
    // e.g.:
    // returnFile: fileActions.retrieveFile
  };

  let workers = {
    upload: createWorker( require.resolve('../jsx/workers/uploader'))
  }

  let middleWares = [
    asyncRouter(actions),
    workRouter(workers)
  ];

  if(process.env.NODE_ENV != 'production') middleWares.push(ReduxLogger.createLogger());

  let store = Redux.applyMiddleware(...middleWares)(
    Redux.createStore
  )(
    Redux.combineReducers(reducers)
  );

  return store;
}


export default createStore;
