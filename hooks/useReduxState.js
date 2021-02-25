import {useContext, useState, useEffect} from 'react';
import {ReactReduxContext} from "react-redux";

const ident = (state) => state;

export function useReduxState(selectorFunction = ident) {
  const { store } = useContext(ReactReduxContext);
  const [state, setState] = useState(() => selectorFunction(store.getState()));

  useEffect(() => {
     return store.subscribe(() => {
      const newState = store.getState();
      const selected = selectorFunction(newState)
      if (!Array.isArray(selected) && typeof selected === 'object') {
        if (!isShallowDifferent(selected, state)) {
          return;
        }
      } else {
        if (selected === state) return;
      }

      setState(selected);
    });
  }, [store]);

  return state;
}

function isShallowDifferent(o1, o2) {
  for (let key in o1) {
    if (o1[key] !== o2[key]) {
      return true;
    }
  }

  for (let key in o2) {
    if (o1[key] !== o2[key]) {
      return true;
    }
  }

  return false;
}

