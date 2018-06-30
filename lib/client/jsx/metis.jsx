import React from 'react';
import ReactDOM from 'react-dom';
import { Provider } from 'react-redux';
import Cookies from 'js-cookie';

import createStore from './store';
import MetisUI from './components/metis-ui';

class Metis {
  constructor() {
    this.store = createStore();

    // add the user
    this.store.dispatch({
      type: 'ADD_USER',
      token: Cookies.get(CONFIG.token_name)
    });

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
