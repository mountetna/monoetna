// a middleware that takes routes matching
// action names and sends data to them.
import { camelCase } from '../utils/format';

const asyncDispatcher = actions => {
  return store => next => action => {
    let { type, ...args } = action;

    let name = camelCase(type);

    // execute the action if it is one of our
    // defined action handlers
    if (actions[name]) {
      actions[name](args)(store.dispatch, store.getState);
      return;
    }

    return next(action);
  }
}

export default asyncDispatcher;
