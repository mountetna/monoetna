// Framework libraries.
import React, {useState, useCallback, useEffect, useMemo} from 'react';
import 'regenerator-runtime/runtime';
import useAsyncWork from 'etna-js/hooks/useAsyncWork';

// Module imports.
import {setLocation} from 'etna-js/actions/location_actions';
import {useReduxState} from 'etna-js/hooks/useReduxState';
import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';

const loadingDiv = (
  <div className='browser'>
    <div id='loader-container'>
      <div className='loader'>Loading...</div>
    </div>
  </div>
);

const errorDiv = (
  <div className='browser'>
    <div id='loader-container'>
      <div className='loader'>Failed to load.</div>
    </div>
  </div>
);

export default function Browser() {
  const invoke = useActionInvoker();
  const browserState = useReduxState(browserStateOf());

  return <div>Stub for Vulcan!</div>;
}

function browserStateOf() {
  return (state) => {
    return {};
  };
}
