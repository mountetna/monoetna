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

  removeBucket() {
    let { removeBucket, bucket } = this.props;

    removeFolder(bucket);
  }

  render() {
    let { bucket } = this.props;
    let items = [
      { label: 'Configure bucket', callback: this.configureBucket.bind(this), show: true },
      { label: 'Remove bucket', callback: this.removeBucket.bind(this) },
    ];
    return <MenuControl items={items}/>;
  }
}

export default connect(
  null,
  (dispatch) => ({
    removeBucket: (bucket) => dispatch({ type: 'REMOVE_BUCKET', bucket }),
    updateBucket: (bucket) => dispatch({ type: 'UPDATE_BUCKET', bucket}),
    showDialog: (dialog) => dispatch({ type: 'SHOW_DIALOG', dialog})
  })
)(BucketControl);

