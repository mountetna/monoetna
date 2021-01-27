import {combineReducers, createStore, applyMiddleware} from 'redux';

import thunk from 'redux-thunk';
import {createLogger} from 'redux-logger';

// Reducers.
import messages from 'etna-js/reducers/message_reducer';
import exchanges from 'etna-js/reducers/exchanges_reducer';
import location from 'etna-js/reducers/location_reducer';
import janus from 'etna-js/reducers/janus-reducer';

export const VulcanStore = () => {
  let reducers = combineReducers({
    messages,
    exchanges,
    location,
    janus
  });

  let middlewares = [thunk];

  if (process.env.NODE_ENV != 'production')
    middlewares.push(createLogger({collapsed: true}));

  return createStore(reducers, {}, applyMiddleware(...middlewares));
};
