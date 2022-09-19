import * as Redux from 'redux';
import thunk from 'redux-thunk';
import * as ReduxLogger from 'redux-logger';
import directory from './reducers/directory-reducer';
import user from './reducers/user-reducer';
import dialog from './reducers/dialog-reducer';
import janus from 'etna-js/reducers/janus-reducer';

import * as fileActions from './actions/file_actions';
import * as folderActions from './actions/folder_actions';
import * as bucketActions from './actions/bucket_actions';
import * as uploadActions from 'etna-js/upload/actions/upload_actions';

import asyncDispatcher from 'etna-js/dispatchers/async-dispatcher';
import workDispatcher from 'etna-js/dispatchers/work-dispatcher';

const createStore = () => {
  let reducers = {
    directory,
    dialog,
    user,
    janus
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

  let middleWares = [thunk, asyncDispatcher(actions), workDispatcher()];

  if (process.env.NODE_ENV != 'production') {middleWares.unshift(ReduxLogger.createLogger({collapsed: true}));}

  return Redux.createStore(
    Redux.combineReducers(reducers),
    {},
    Redux.applyMiddleware(...middleWares)
  );
};

export default createStore;
