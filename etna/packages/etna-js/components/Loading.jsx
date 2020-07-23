import React, { useState, useEffect } from 'react';

/*
  Convenience component that will conditionally replace the contents with the 'loading' value provided when it is
  truthy, optionally after a given delay count in ms.
  Alternatively, if loading is simply the boolean state true, it will simply return null for convenience.

  Usage generally involves passing the react fragment into loading like so:
  <Loading loading={loading && <MyLoadingBar/>} delay={200}>
    Not loading content here!
  </Loading>
 */
export function Loading({ children, loading = null, delay = 0, cacheLastView = false }) {
  const [hasDelayed, setHasDelayed] = useState(true);
  const [lastChildren, setLastChildren] = useState(children);

  useEffect(() => {
    if (!loading) setLastChildren(children);
    if (!loading || !delay) return;
    setHasDelayed(false);
    const timer = setTimeout(() => setHasDelayed(true), delay);
    return () => clearTimeout(timer);
  }, [ loading, setHasDelayed, delay ]);

  if (loading && hasDelayed) {
    if (loading === true) return null;
    return loading;
  }

  if (loading && !hasDelayed && cacheLastView) {
    return lastChildren;
  }

  return <React.Fragment>
    {children}
  </React.Fragment>;
}