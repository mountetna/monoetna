import * as Redux from 'redux';
import * as ReduxLogger from 'redux-logger';
import directory from './reducers/directory-reducer';
import user from './reducers/user-reducer';
import dialog from './reducers/dialog-reducer';

import * as fileActions from './actions/file_actions';
import * as folderActions from './actions/folder_actions';
import * as bucketActions from './actions/bucket_actions';
import * as uploadActions from './actions/upload_actions';

import asyncRouter from './routers/async-router';
import workRouter from './routers/work-router';

const createStore = () => {
  let reducers = {
    directory,
    dialog,
    user
  };

  // action handlers to import
  let actions = {
    ...fileActions,
    ...folderActions,
    ...bucketActions,
    ...uploadActions

    // here you may define aliases to other actions,
    // e.g.:
    // returnFile: fileActions.retrieveFile
  };


  let middleWares = [
    asyncRouter(actions),
    workRouter()
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
