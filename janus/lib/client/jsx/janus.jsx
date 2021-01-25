import React from 'react';
import ReactDOM from 'react-dom';
import { Provider } from 'react-redux';
import Cookies from 'js-cookie';

import createStore from './store';
import JanusUI from './janus-ui';

function Janus() {
  this.store = createStore();

  // add the user
  this.store.dispatch({
    type: 'ADD_USER',
    token: Cookies.get(CONFIG.token_name)
  });

  // build the UI
  ReactDOM.render(
    <Provider store={ this.store }>
      <JanusUI/>
    </Provider>,
    document.getElementById('root')
  );
}

window.janus = new Janus();
