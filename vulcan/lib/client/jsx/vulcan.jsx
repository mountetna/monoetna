// Framework Libraries.
import 'regenerator-runtime/runtime';
import * as React from 'react';
import * as ReactDOM from 'react-dom';

import {Provider} from 'react-redux';
import VulcanUI from './vulcan_ui';
import {VulcanStore} from './vulcan_store';
import * as Cookies from 'etna-js/utils/cookies';

class VulcanApplication {
  constructor(props, container_id) {
    // create the store
    this.createStore();

    // add user info from the token to the store
    this.store.dispatch({
      type: 'ADD_TOKEN_USER',
      token: Cookies.getItem(CONFIG.token_name)
    });

    // create the base component
    this.createUI(props, container_id);
  }

  createStore() {
    this.store = VulcanStore();
  }

  createUI({environment}, container_id) {
    ReactDOM.render(
      <Provider store={this.store}>
        <VulcanUI
          environment={environment}
          path={decodeURI(window.location.pathname)}
        />
      </Provider>,
      document.getElementById(container_id)
    );
  }
}

window.VulcanApp = VulcanApplication;
