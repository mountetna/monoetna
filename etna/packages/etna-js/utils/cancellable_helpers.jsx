import 'regenerator-runtime';
import {Cancellable, cancelledAsMaybe} from "./cancellable";
import {useState, useEffect, useCallback} from 'react';

export function useAsyncCallback(fn, deps, cleanup = () => null) {
  const context = useContext(new Cancellable());
  const cancel = useCallback(() => {
    context.cancellable.cancel();
    const newCancellable = context.cancellable = new Cancellable();
    return newCancellable;
  }, []);

  useEffect(() => {
    return () => {
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
    context.run(fn()).finally(cleanup).then(cancelledAsMaybe).then(setResult, (e) => setError([e]));
  }, deps);

  return [result, error];
}

// Utility functions leveraged by typescript when handling generators or promises
// inside of generic generators.
export function* runAsync(fn) {
  const result = yield fn();
  return result;
}

export function* runGen(gen) {
  return yield* gen;
}

export function* runPromise(v) {
  const result = yield v;
  return result;
}
