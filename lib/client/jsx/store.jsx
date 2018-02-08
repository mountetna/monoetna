import * as Redux from 'redux';
import * as ReduxLogger from 'redux-logger';
import asyncRouter from './async-router';
import workRouter from './work-router';

const createStore = (reducers, actions, workers) => {
  let reducer = Redux.combineReducers(reducers);

  let middleWares = [
    asyncRouter(actions),
    workRouter(workers)
  ];

  if(process.env.NODE_ENV != 'production') middleWares.push(ReduxLogger.createLogger());

  let store = Redux.applyMiddleware(...middleWares)(Redux.createStore)(reducer);

  return store;
}

export default createStore;
