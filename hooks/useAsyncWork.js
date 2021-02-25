import {useState, useCallback, useEffect} from 'react';
import {Cancellable} from "../utils/cancellable";

/*
  A convenience hook that wraps an async function or cancellable protocol supporting generator function such that
  'loading' state is tracked automatically, and possibly such that 'cancellability' is tracked via watching a list of
  values for changes.

  IFF f is just an async function, returns a pair [loading, wrapper] such that
    loading is true after any invocation of wrapper until the promise returned by that wrapped rejects or resolves
    wrapper opaquely passes its arguments to f and returns its values, inserting hooks into the future promise
    to handle the loading state only.

  IFF f is a generator AND cancelWhenChange is some list, returns a pair [loading, wrapper] such that
    loading is true after any invocation of wrapper until the generator returned by wrapper completes or is cancelled
    wrapper opaquely passes its arguments to f and uses it's generation as an argument to Cancellable.run -- meaning it
    should yield promises to perform asynchronous work.  The resulting promise of wrapper will resolve with the pair
    { result, cancelled } just as Cancellable.run is called.  loading will be managed as a natural side effect of that
    promise ending, either via cancellation or completion.  Cancellation is handled via useEffect -- any further renders
    that result in calls to this hook may trigger a cancellation of the last known invocation of f.  Any consecutive
    invocations of f cancel any previous invocations.
 */
export default function useAsyncWork(f, { cancelWhenChange = undefined, renderedState = undefined }) {
  const [loading, setLoading] = useState(false);
  const [lastPendingResolve, setPendingResolve] = useState(null);
  if (lastPendingResolve) {
    lastPendingResolve(renderedState);
    setPendingResolve(null);
  }
  // const [last, setLoading] = useState(false);

  const [cancellable, setCancellable] = cancelWhenChange ? useState(new Cancellable()) : [null, null];

  if (cancellable) useEffect(function () {
    // Return the cancel as a cleanup operation -- this will get invokved only between the cancelWhenChange has new
    // values, or the component is dismounted.
    return () => cancellable.cancel();
  }, cancelWhenChange);

  const wrapper = useCallback(_wrapper, [setLoading, f, cancellable, setCancellable]);
  const awaitNextRender = useCallback(_awaitNextRender, [setPendingResolve]);
  return [loading, wrapper, awaitNextRender];

  function _awaitNextRender() {
    return new Promise((resolve) => {
      setPendingResolve(() => resolve);
    })
  }

  function _wrapper() {
    let result = f.apply(this, arguments);
    setLoading(true);


    if (cancellable) {
      cancellable.cancel();
      const newCancellable = new Cancellable();
      setCancellable(newCancellable);
      result = newCancellable.run(result);
    }

    return result.finally(() => setLoading(false));
  }
}