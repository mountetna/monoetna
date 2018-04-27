import React from 'react';
import ReactDOM from 'react-dom';
import { Provider } from 'react-redux';
import Cookies from 'js-cookie';

import metisStore from './store';
import MetisUI from './components/metis-ui';

import files from './reducers/files-reducer';
import user from './reducers/janus-log-reducer';

import * as authActions from './actions/auth_actions';
import * as fileActions from './actions/file_actions';
import * as uploadActions from './actions/upload_actions';
import * as userActions from './actions/user_actions';

import { createWorker } from './workers';

const createStore = () => {
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

class Metis {
  constructor() {
    this.store = createStore();

    // add the user
    this.store.dispatch({
      type: 'ADD_USER',
      token: Cookies.get(CONFIG.token_name)
    });

    // request an initial download view
    this.store.dispatch({ type: 'RETRIEVE_FILES' });

    // build the UI
    ReactDOM.render(
      <Provider store={ this.store }>
        <MetisUI />
      </Provider>,
      document.getElementById('ui-group')
    );
  }
}

let metis = new Metis();
