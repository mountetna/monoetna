import * as React from 'react';
import { connect } from 'react-redux';

import ListBody  from './list/list-body';
import ListHead  from './list/list-head';
import FolderBreadcrumb from './folder-breadcrumb';
import ControlBar from './control-bar';

const COLUMNS = [
  { name: 'type', width: '90px' },
  { name: 'name', width: '60%' },
  { name: 'status', width: '90px', hide: true },
  { name: 'updated', width: '30%' },
  { name: 'size', width: '10%' },
  { name: 'control', width: '100px', hide: true }
];

const COLUMN_WIDTHS = COLUMNS.reduce( (widths,column) => {
  widths[column.name] = column.width;
  return widths;
}, {} );

const INVALID = '\ninvalid\n';

const InvalidFolder = () => <div className='invalid-folder-view-group'>
  Invalid folder!
</div>;

class FolderView extends React.Component {
  componentDidMount() {
    let { bucket_name, folder_name, retrieveFiles } = this.props;

    retrieveFiles(bucket_name, folder_name);
  }

  render() {
    let { bucket_name, folder_name } = this.props;
    if (folder_name == INVALID) return <InvalidFolder/>

    return (
      <div className='folder-view-group'>
        <div className='control-group'>
          <FolderBreadcrumb folder_name={folder_name} bucket_name={bucket_name}/>
          <ControlBar bucket_name={ bucket_name } folder_name={ folder_name }/>
        </div>
        <div className='listing-group'>
          <ListHead columns={ COLUMNS }/>
          <ListBody widths={ COLUMN_WIDTHS } folder_name={folder_name} bucket_name={bucket_name}/>
        </div>
      </div>
    );
  }
}

export default connect(
  // map state
  null,

  // map dispatch
  (dispatch) => ({
    retrieveFiles: (bucket_name, folder_name) => dispatch({type: 'RETRIEVE_FILES', bucket_name, folder_name})
  })
)(FolderView);
