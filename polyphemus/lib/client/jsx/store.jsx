import * as Redux from 'redux';
import thunk from 'redux-thunk';
import * as ReduxLogger from 'redux-logger';
import user from 'etna-js/reducers/user-reducer';

const createStore = () => {
  let reducers = { user };

  let middleWares = [thunk];

  if (process.env.NODE_ENV != 'production')
    middleWares.unshift(ReduxLogger.createLogger({collapsed: true}));

  return Redux.createStore(
    Redux.combineReducers(reducers),
    {},
    Redux.applyMiddleware(...middleWares)
  );
};

export default createStore;
