import * as Redux from 'redux';
import thunk from 'redux-thunk';
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
  const reducers = {
    user,
    janus,
    location,
    rules: rulesReducer,
    names: namesReducer,
  };

  const middleWares = [thunk];

  let composeEnhancers = Redux.compose;

  // @ts-ignore
  if (process.env.NODE_ENV != 'production') {
    // middleWares.unshift(ReduxLogger.createLogger({ collapsed: true }));

    // @ts-ignore
    if (typeof window === 'object' && window.__REDUX_DEVTOOLS_EXTENSION_COMPOSE__) {
      // @ts-ignore
      composeEnhancers = window.__REDUX_DEVTOOLS_EXTENSION_COMPOSE__({
        // options: https://github.com/zalmoxisus/redux-devtools-extension/blob/master/docs/API/Arguments.md
        // trace: true,
      });
    }
  }

  const enhancer = composeEnhancers(
    Redux.applyMiddleware(...middleWares),
    // other store enhancers if any
  );

  return Redux.createStore(
    Redux.combineReducers(reducers),
    {},
    enhancer,
  );
};

export default createStore;
