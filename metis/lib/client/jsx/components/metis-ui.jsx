import React from 'react';

import MetisNav from './metis-nav';
import RootView from 'etna-js/components/RootView';
import FolderView from './folder-view';
import BucketView from './bucket-view';
import ModalDialog from './modal-dialog';

import {findRoute, setRoutes} from 'etna-js/dispatchers/router';

import {ThemeProvider} from '@material-ui/core/styles';
import {createEtnaTheme} from 'etna-js/style/theme';

import {IngestProvider} from '../contexts/ingest_context';

const theme = createEtnaTheme('goldenrod', '#066306');

const ROUTES = [
  {
    name: 'root',
    template: '',
    component: RootView
  },
  {
    name: 'bucket',
    template: ':project_name',
    component: BucketView
  },
  {
    name: 'folder',
    template: ':project_name/browse/:bucket_name/*folder_name',
    component: FolderView
  },
  {
    name: 'root_folder',
    template: ':project_name/browse/:bucket_name',
    component: FolderView
  }
];

setRoutes(ROUTES);

const Invalid = () => <div>Path invalid</div>;

const MetisUI = () => {
  let {route, params} = findRoute({path: window.location.pathname}, ROUTES);
  let Component = route ? route.component : Invalid;

  return (
    <div id='metis-group'>
      <IngestProvider>
        <ThemeProvider theme={theme}>
          <MetisNav />
          <Component {...params} />
          <ModalDialog />
        </ThemeProvider>
      </IngestProvider>
    </div>
  );
};

export default MetisUI;
