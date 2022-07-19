import React from 'react';
import ReactDOM from 'react-dom';
import { Provider } from 'react-redux';
import Cookies from 'js-cookie';

import createStore from './store';
import PolyphemusUI from './polyphemus-ui';

function Polyphemus() {
  const store = createStore();

  // add the user
  store.dispatch({
    type: 'ADD_USER',
    token: Cookies.get(CONFIG.token_name)
  });

  // build the UI
  ReactDOM.render(
    <Provider store={ store }>
      <PolyphemusUI/>
    </Provider>,
    document.getElementById('root')
  );
}

window.polyphemus = new Polyphemus();
