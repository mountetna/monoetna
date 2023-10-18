import * as Redux from 'redux';
import thunk from 'redux-thunk';
import * as ReduxLogger from 'redux-logger';
import user from 'etna-js/reducers/user-reducer';
import janus from 'etna-js/reducers/janus-reducer';
import location from 'etna-js/reducers/location_reducer';
import { RulesState, rulesReducer } from './reducers/rules';
import { NamesState, namesReducer } from './reducers/names';



export interface LocationState {
  path: string
  search: URLSearchParams
  hash: string | null
}


export interface State {
  user: Record<string, any>
  janus: { projects: Record<any, any>[] }
  location: LocationState
  rules: RulesState
  names: NamesState
}


const createStore = () => {
  let reducers = {
    user,
    janus,
    location,
    rules: rulesReducer,
    names: namesReducer,
  };

  let middleWares = [thunk];

  if (process.env.NODE_ENV != 'production') {
    middleWares.unshift(ReduxLogger.createLogger({ collapsed: true }));
  }

  return Redux.createStore(
    Redux.combineReducers(reducers),
    {},
    Redux.applyMiddleware(...middleWares)
  );
};

export default createStore;
