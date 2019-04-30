import * as React from 'react';
import { connect } from 'react-redux';

import MetisNav from './metis-nav';
import FolderView from './folder-view';
import BucketView from './bucket-view';
import ModalDialog from './modal-dialog';

import { findRoute, setRoutes } from '../router';

const ROUTES = [
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

class MetisUI extends React.Component {
  render() {
    let { route, params }  = findRoute({ path: window.location.pathname } ,ROUTES);
    let Component = route ? route.component : Invalid;

    return (
      <div id='metis-group'>
        <MetisNav/>
        <Component {...params}/>
        <ModalDialog/>
      </div>
    );
  }
}

export default MetisUI;
