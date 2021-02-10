// Framework libraries.
import React, {useState, useCallback, useEffect, useMemo} from 'react';
import 'regenerator-runtime/runtime';
import useAsyncWork from 'etna-js/hooks/useAsyncWork';

// Module imports.
import {setLocation} from 'etna-js/actions/location_actions';
import {useReduxState} from 'etna-js/hooks/useReduxState';
import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';

import {selectProjects} from 'etna-js/selectors/janus-selector';
import {projectNameFull} from 'etna-js/utils/janus';
import {fetchProjectsAction} from 'etna-js/actions/janus-actions';

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

  let {projects} = browserState;

  useEffect(() => {
    invoke(fetchProjectsAction());
  });

  return (
    <div>
      Select a workflow for{' '}
      {projectNameFull(projects, CONFIG.project_name) || CONFIG.project_name}
    </div>
  );
}

function browserStateOf() {
  return (state) => {
    return {projects: selectProjects(state)};
  };
}
