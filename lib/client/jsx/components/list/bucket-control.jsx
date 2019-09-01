import * as React from 'react';
import { connect } from 'react-redux';
import { copyText } from '../../utils/copy';
import { filePath } from '../../utils/file';
import MenuControl from '../menu-control';

class BucketControl extends React.Component {
  configureBucket() {
    let { showDialog, updateBucket, bucket } = this.props;
    let { bucket_name, description, access } = bucket;

    console.log(bucket);
    let dialog = {
      type: 'configure-bucket',
      width: 300,
      height: 1,
      top: 0, right: 0, left: 0, bottom: 0,
      updateBucket,
      bucket_name, description, access
    }
    showDialog(dialog);
  }

  destroyBucket() {
    let { destroyBucket, bucket } = this.props;

    destroyBucket(bucket);
  }

  render() {
    let { bucket } = this.props;
    let items = [
      { label: 'Configure bucket', callback: this.configureBucket.bind(this), show: true },
      bucket.count == 0 && { label: 'Remove bucket', callback: this.destroyBucket.bind(this) }
    ].filter(_=>_);
    return <MenuControl items={items}/>;
  }
}

export default connect(
  null,
  (dispatch) => ({
    destroyBucket: (bucket) => dispatch({ type: 'DESTROY_BUCKET', bucket }),
    updateBucket: (bucket) => dispatch({ type: 'UPDATE_BUCKET', bucket}),
    showDialog: (dialog) => dispatch({ type: 'SHOW_DIALOG', dialog})
  })
)(BucketControl);

