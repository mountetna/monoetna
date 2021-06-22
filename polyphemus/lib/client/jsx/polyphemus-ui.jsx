import * as React from 'react';
import { ThemeProvider } from '@material-ui/core/styles';

import PolyphemusNav from './polyphemus-nav';
import PolyphemusMain from './polyphemus-main';

import { findRoute, setRoutes } from 'etna-js/dispatchers/router';

import { createEtnaTheme } from 'etna-js/style/theme';

const theme = createEtnaTheme("#3684fd","#77c");

const ROUTES = [
  {
    name: 'main',
    template: '',
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
