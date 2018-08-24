import * as React from 'react';
import { connect } from 'react-redux';
import { copyText } from '../../utils/copy';
import { filePath } from '../../utils/file';

class FileControl extends React.Component{
  constructor(props) {
    super(props);
    this.state = { dialog: false };
  }

  unprotectFile() {
    let { unprotectFile, file} = this.props;

    unprotectFile(file);
  }

  protectFile() {
    let { protectFile, file} = this.props;

    protectFile(file);
  }

  renameFile() {
    let { renameFile, current_folder, file } = this.props;
    let new_file_name = prompt("What is the new name of this file?", file.file_name);
    if (new_file_name) renameFile(file, filePath(current_folder, new_file_name));
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

    let items = file.read_only ?
      [
        { label: 'Unprotect file', callback: this.unprotectFile.bind(this) }
      ] : [
        { label: 'Rename file', callback: this.renameFile.bind(this) },
        { label: 'Protect file', callback: this.protectFile.bind(this) },
        { label: 'Remove file', callback: this.removeFile.bind(this) },
      ];

    let dialog = {
      type: 'list',
      items: [
        { label: 'Copy download link', callback: this.copyLink.bind(this) },
        ...items
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
    unprotectFile: (file) => dispatch({ type: 'UNPROTECT_FILE', file }),
    protectFile: (file) => dispatch({ type: 'PROTECT_FILE', file }),
    renameFile: (file, new_file_path) => dispatch({ type: 'RENAME_FILE', file, new_file_path })
  })
)(FileControl);

