import * as React from 'react';
import { ThemeProvider } from '@material-ui/core/styles';

import GnomonNav from './gnomon-nav';
import GnomonMain from './gnomon-main';
import DecomposeIdentifier from './decompose-identifier';
import ComposeIdentifier from './compose-identifier';
import RuleEditor from './rule-editor';
import RootView from 'etna-js/components/RootView';

import { findRoute, setRoutes } from 'etna-js/dispatchers/router';

import { createEtnaTheme } from 'etna-js/style/theme';

const theme = createEtnaTheme("#444","#888");

const ROUTES = [
  {
    template: '',
    component: RootView
  },
  {
    template: ':project_name/',
    component: GnomonMain
  },
  {
    template: ':project_name/identify/*identifier',
    component: DecomposeIdentifier
  },
  {
    template: ':project_name/create/:rule_name',
    component: ComposeIdentifier
  },
  {
    template: ':project_name/rules',
    component: RuleEditor
  }
];

setRoutes(ROUTES);

const Invalid = () => <div>Path invalid</div>;

const GnomonUI = () => {
  let { route, params }  = findRoute({ path: window.location.pathname }, ROUTES);
  let Component = route ? route.component : Invalid;

  return (
    <ThemeProvider theme={theme}>
      <div id='gnomon-group'>
        <GnomonNav/>
        <Component {...params}/>
      </div>
    </ThemeProvider>
  );
}

export default GnomonUI;
