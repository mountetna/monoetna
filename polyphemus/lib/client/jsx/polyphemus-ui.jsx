import * as React from 'react';
import { ThemeProvider } from '@material-ui/core/styles';

import PolyphemusNav from './polyphemus-nav';
import PolyphemusMain from './polyphemus-main';
import RootView from 'etna-js/components/RootView';

import { findRoute, setRoutes } from 'etna-js/dispatchers/router';

import { createEtnaTheme } from 'etna-js/style/theme';

const theme = createEtnaTheme("#a8ed44","#d18e47");

const ROUTES = [
  {
    template: '',
    component: RootView
  },
  {
    template: ':project_name/',
    component: PolyphemusMain
  }
];

setRoutes(ROUTES);

const Invalid = () => <div>Path invalid</div>;

const PolyphemusUI = () => {
  let { route, params }  = findRoute({ path: window.location.pathname }, ROUTES);
  let Component = route ? route.component : Invalid;

  return (
    <ThemeProvider theme={theme}>
      <div id='polyphemus-group'>
        <PolyphemusNav/>
        <Component {...params}/>
      </div>
    </ThemeProvider>
  );
}

export default PolyphemusUI;
