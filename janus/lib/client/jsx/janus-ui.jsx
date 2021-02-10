import * as React from 'react';

import JanusNav from './janus-nav';
import JanusMain from './janus-main';
import JanusAdmin from './janus-admin';
import ProjectView from './project-view';

import { findRoute, setRoutes } from 'etna-js/dispatchers/router';

const ROUTES = [
  {
    name: 'main',
    template: '',
    component: JanusMain
  },
  {
    name: 'admin',
    template: 'admin',
    component: JanusAdmin
  },
  {
    name: 'projects',
    template: ':project_name',
    component: ProjectView
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
