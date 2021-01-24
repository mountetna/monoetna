import * as React from 'react';

import JanusNav from './janus-nav';
import JanusRoot from './janus-root';

import { findRoute, setRoutes } from 'etna-js/dispatchers/router';

const ROUTES = [
  {
    name: 'root',
    template: '',
    component: JanusRoot
  }
];

setRoutes(ROUTES);

const Invalid = () => <div>Path invalid</div>;

const JanusUI = () => {
  let { route, params }  = findRoute({ path: window.location.pathname }, ROUTES);
  let Component = route ? route.component : Invalid;

  return (
    <div id='janus-group'>
      <JanusNav/>
      <Component {...params}/>
    </div>
  );
}

export default JanusUI;
