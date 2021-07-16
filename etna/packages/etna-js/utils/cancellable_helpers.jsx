import 'regenerator-runtime';
import {Cancellable, cancelledAsMaybe} from "./cancellable";
import {useState, useEffect, useCallback} from 'react';

export function useAsyncCallback(fn, deps, cleanup = () => null) {
  const [_, setContext] = useState(() => new Cancellable());
  const cancel = useCallback(() => {
    const newCancellable = new Cancellable();
    setContext(running => {
      running.cancel();
      return newCancellable;
    })

    return newCancellable;
  }, []);

  useEffect(() => {
    return () => {
      console.log('cancelling from efect')
      cancel();
    }
  }, deps);

  const start = useCallback(function (...args) {
    const next = cancel();
    return next.run(fn(...args)).finally(cleanup);
  }, deps);

  return [start, cancel];
}

export function useWithContext(fn, deps) {
  useEffect(() => {
    const cancellable = new Cancellable();
    fn(cancellable);
    return () => cancellable.cancel();
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, deps);
}

export function useAsync(fn, deps= [], cleanup = () => null) {
  const [result, setResult] = useState(null);
  const [error, setError] = useState(null);

  useWithContext((context) => {
    setError(null);
    setResult(null);
    context.run(fn()).then(cancelledAsMaybe).then(setResult, (e) => setError([e])).finally(cleanup);
  }, deps);

  return [result, error];
}

// Utility functions leveraged by typescript
export function* runAsync(fn) {
  const result = yield fn();
  return result;
}

export function* runPromise(v) {
  const result = yield v;
  return result;
}
