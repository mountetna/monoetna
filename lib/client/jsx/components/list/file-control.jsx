import * as React from 'react';
import { connect } from 'react-redux';
import { copyText } from '../../utils/copy';

class FileControl extends React.Component{
  constructor(props) {
    super(props);
    this.state = { dialog: false };
  }

  protectFile() {
    let { protectFile, file} = this.props;

    protectFile(file);
  }


  renameFile() {
    alert(`renaming file ${this.props.file.file_name}`);
  }

  removeFile() {
    let { removeFile, file} = this.props;

    removeFile(file);
  }

  copyLink() {
    let { file: { download_url } } = this.props;

    copyText(download_url);
  }

  dismissDialog() {
    this.setState({dialog: false});
  }

  showDialog() {
    let { showDialog, file } = this.props;
    let { top, left } = this.control.getBoundingClientRect();

    this.setState({ dialog: true });

    let dialog = {
      type: 'list',
      items: [
        { label: 'Copy download link', callback: this.copyLink.bind(this) },
        { label: 'Rename file', callback: this.renameFile.bind(this) },
        { label: 'Protect file', callback: this.protectFile.bind(this) },
        { label: 'Remove file', callback: this.removeFile.bind(this) },
      ],
      dismiss: this.dismissDialog.bind(this),
      top,
      left
    };

    showDialog(dialog);
  }

  setControlRef(ref) {
    this.control = ref;
  }

  render() {
    let { dialog } = this.state;
    let className = `control-btn-group ${dialog ? 'selected' : ''}`;
    return (
      <div ref={ this.setControlRef.bind(this) } className={className} onClick={this.showDialog.bind(this)}>
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
    showDialog: (dialog) => dispatch({ type: 'SHOW_DIALOG', dialog}),
    removeFile: (file) => dispatch({ type: 'REMOVE_FILE', file }),
    protectFile: (file) => dispatch({ type: 'PROTECT_FILE', file })
  })
)(FileControl);

