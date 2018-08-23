import * as React from 'react';
import { connect } from 'react-redux';

class FileControl extends React.Component{
  protectFile() {
  }

  showDialog() {
    let { showDialog, file } = this.props;
    let { top, left } = this.control.getBoundingClientRect();

    let dialog = {
      type: 'list',
      items: [
        { label: 'Copy download link', callback: this.protectFile.bind(this) },
        { label: 'Rename file', callback: this.protectFile.bind(this) },
        { label: 'Protect file', callback: this.protectFile.bind(this) },
        { label: 'Remove file', callback: this.protectFile.bind(this) },
      ],
      top,
      left
    };

    showDialog(dialog);
  }

  setControlRef(ref) {
    this.control = ref;
  }

  render() {
    return (
      <div ref={ this.setControlRef.bind(this) } className='file-control-group' onClick={this.showDialog.bind(this)}>
        &bull;
        &bull;
        &bull;
      </div>
    )
  }
}

export default connect(
  null,
  (dispatch) => ({
    showDialog: (dialog) => dispatch({ type: 'SHOW_DIALOG', dialog})
  })
)(FileControl);

