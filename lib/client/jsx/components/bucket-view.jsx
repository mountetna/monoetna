import * as React from 'react';
import { connect } from 'react-redux';

import ListBucket from './list/list-bucket';
import ListHead  from './list/list-head';
import FolderBreadcrumb from './folder-breadcrumb';
import { ControlButton } from './control-bar';

import { selectBuckets } from '../selectors/directory-selector';

const COLUMNS = [
  { name: 'type', width: '90px' },
  { name: 'name', width: '60%' },
  { name: 'description', width: '60%' },
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

class BucketView extends React.Component {
  componentDidMount() {
    let { bucket_name, folder_name, retrieveBuckets } = this.props;

    retrieveBuckets();
  }

  render() {
    let { project_name, buckets, showDialog, createBucket } = this.props;

    let dialog = {
      type: 'create-bucket',
      width: 300,
      height: 1,
      top: 0, right: 0, left: 0, bottom: 0,
      createBucket
    }

    return (
      <div className='bucket-view-group'>
        <div className='control-group'>
          <FolderBreadcrumb/>
          <div id='control-bar'>
            <ControlButton
              onClick={ () => showDialog(dialog) }
              title='Create bucket'
              icon='trash'
              overlay='plus'/>
          </div>
        </div>
        <ListHead columns={ COLUMNS }/>
        <div id='list-body-group'>
        {
          Object.values(buckets).sort(b=>b.created_at).map(
            bucket => <ListBucket key={bucket.bucket_name} widths={COLUMN_WIDTHS} bucket={bucket}/>
          )
        }
        </div>
      </div>
    );
  }
}

export default connect(
  // map state
  (state) => ({
    buckets: selectBuckets(state)
  }),
  // map dispatch
  (dispatch) => ({
    retrieveBuckets: (bucket_name, folder_name) => dispatch({type: 'RETRIEVE_BUCKETS' }),
    createBucket: (bucket) => dispatch({type: 'CREATE_BUCKET', bucket }),
    showDialog: (dialog) => dispatch({ type: 'SHOW_DIALOG', dialog})
  })
)(BucketView);
