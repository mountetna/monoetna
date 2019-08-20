import * as React from 'react';
import { connect } from 'react-redux';

import ListBucket from './list/list-bucket';
import ListHead  from './list/list-head';
import FolderBreadcrumb from './folder-breadcrumb';

import { selectBuckets } from '../selectors/directory-selector';

const COLUMNS = [
  { name: 'type', width: '90px' },
  { name: 'name', width: '60%' },
  { name: 'description', width: '40%' }
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
    let { project_name } = this.props;
    let { buckets } = this.props;

    return (
      <div className='bucket-view-group'>
        <div className='control-group'>
          <FolderBreadcrumb/>
        </div>
        <ListHead columns={ COLUMNS }/>
        <div id='list-body-group'>
        {
          buckets.map(
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
    retrieveBuckets: (bucket_name, folder_name) => dispatch({type: 'RETRIEVE_BUCKETS' })
  })
)(BucketView);
