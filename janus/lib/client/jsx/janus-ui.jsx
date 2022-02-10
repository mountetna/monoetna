import * as React from 'react';
import { ThemeProvider } from '@material-ui/core/styles';

import JanusNav from './janus-nav';
import JanusMain from './janus-main';
import JanusAdmin from './janus-admin';
import JanusSettings from './janus-settings';
import ProjectView from './project-view';
import FlagsView from './flags/flags-view';

import { findRoute, setRoutes } from 'etna-js/dispatchers/router';

import {Notifications} from 'etna-js/components/Notifications';
import Messages from 'etna-js/components/messages';
import { createEtnaTheme } from 'etna-js/style/theme';

const theme = createEtnaTheme("#3684fd","#77c");

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
    name: 'settings',
    template: 'settings',
    component: JanusSettings
  },
  {
    name: 'flags',
    template: 'flags',
    component: FlagsView
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
    <ThemeProvider theme={theme}>
      <div id='janus-group'>
        <Notifications />
        <JanusNav/>
        <Messages />
        <Component {...params}/>
      </div>
    </ThemeProvider>
  );
}

export default JanusUI;
