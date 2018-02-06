import * as Redux from 'redux';
import * as ReduxLogger from 'redux-logger';
import asyncRouter from './async-router';

const createStore = (reducers, actions) => {
  let reducer = Redux.combineReducers(reducers);

  let middleWares = [
    asyncRouter(actions)
  ];

  if(process.env.NODE_ENV != 'production') middleWares.push(ReduxLogger.createLogger());

  let store = Redux.applyMiddleware(...middleWares)(Redux.createStore)(reducer);

  return store;
}

export default createStore;
