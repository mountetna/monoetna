import {useCallback, useContext} from 'react';
import {ReactReduxContext} from "react-redux";

// Allows to avoid having to use the redux-react `connect`, more easily directly invoking actions.
// Supports function (thunk) actions and raw payloads.
// Usage:
// import { doMyAction } from './actions';
// function MyComponent({ blah, moo, ...props }) {
//   const invoke = useActionInvoker();
//   const callback = () => invoke(doMyAction(blah, moo));
// }
export function useActionInvoker() {
  const { store } = useContext(ReactReduxContext);
  return useCallback((action) => {
    if (typeof action === 'function') {
      return action(store.dispatch);
    }

    return store.dispatch(action);
  }, [store]);
}
