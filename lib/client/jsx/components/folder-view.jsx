import * as React from 'react';

import ListBody  from './list/list-body';
import ListHead  from './list/list-head';
import FolderBreadcrumb from './folder-breadcrumb';
import ControlBar from './control-bar';


const COLUMN_WIDTHS = {
  type: '90px',
  name: '60%',
  status: '90px',
  updated: '30%',
  size: '10%',
  control: '100px'
};

const INVALID = '\ninvalid\n';

const FolderView = ({folder}) =>
  folder == INVALID ?
    <div className='invalid-folder-view-group'>
      Invalid folder!
    </div>
  :
    <div className='folder-view-group'>
      <div className='control-group'>
        <FolderBreadcrumb/>
        <ControlBar/>
      </div>
      <div className='listing-group'>
        <ListHead widths={ COLUMN_WIDTHS } />
        <ListBody widths={ COLUMN_WIDTHS }/>
      </div>
    </div>;

export default FolderView;
