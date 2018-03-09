import React from 'react';
import ReactDOM from 'react-dom';
import { Provider } from 'react-redux';
import Cookies from 'js-cookie';

import metisStore from './store';
import MetisUIContainer from './components/metis-ui-container';

import files from './reducers/files-reducer';
import user from './reducers/janus-log-reducer';

import * as authActions from './actions/auth_actions';
import * as fileActions from './actions/file_actions';
import * as uploadActions from './actions/upload_actions';
import * as userActions from './actions/user_actions';

import { createWorker } from './workers';

class MetisUploader {
  constructor() {
    this.store = this.createStore();

    this.updateUser();
    this.buildUI();
    this.retrieveFiles();
  }

  updateUser() {
    this.store.dispatch({
      type: 'ADD_USER',
      token: Cookies.get(CONFIG.token_name)
    });
  }

  createStore() {
    // these are (state, action) => new_state
    let reducers = {
      files,
      user
    };

    // action handlers to import
    let actions = {
      ...authActions,
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

    return metisStore(reducers, actions, workers);
  }

  buildUI(){
    ReactDOM.render(
      <Provider store={ this.store }>
        <MetisUIContainer />
      </Provider>,
      document.getElementById('ui-group')
    );
  }

  retrieveFiles(){
    this.store.dispatch({ type: 'RETRIEVE_FILES' });
  }
}

// Initilize the class.
let metisUploader = new MetisUploader();
