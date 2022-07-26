import React from 'react';
import ReactDOM from 'react-dom';
import { Provider } from 'react-redux';
import Cookies from 'js-cookie';

import createStore from './store';
import GnomonUI from './gnomon-ui';

function Gnomon() {
  const store = createStore();

  // add the user
  store.dispatch({
    type: 'ADD_USER',
    token: Cookies.get(CONFIG.token_name)
  });

  // build the UI
  ReactDOM.render(
    <Provider store={ store }>
      <GnomonUI/>
    </Provider>,
    document.getElementById('root')
  );
}

window.gnomon = new Gnomon();
