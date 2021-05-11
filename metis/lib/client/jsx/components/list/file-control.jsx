import * as React from 'react';
import { connect } from 'react-redux';
import { copyText } from 'etna-js/utils/copy';
import { filePath } from 'etna-js/utils/file';
import MenuControl from '../menu-control';

class FileControl extends React.Component{
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

  copyMetisPath() {
    let { file: { project_name, bucket_name, file_path } } = this.props;

    copyText(`metis://${project_name}/${bucket_name}/${file_path}`);
  }

  downloadFile() {
    let { file: { file_name, download_url } } = this.props;

    let download = document.createElement('a');
    download.setAttribute('href', download_url);
    download.setAttribute('download', file_name);

    download.style.display = 'none';
    document.body.appendChild(download);

    download.click();

    document.body.removeChild(download);
  }

  render() {
    let { file } = this.props;
    let items = [
      { label: 'Download file', callback: this.downloadFile.bind(this) },
      { label: 'Copy download link', callback: this.copyLink.bind(this) },
      { label: 'Copy metis path', callback: this.copyMetisPath.bind(this) }
    ].concat(
      file.read_only ?
      [
        { label: 'Unprotect file', callback: this.unprotectFile.bind(this), role: 'administrator' }
      ] : [
        { label: 'Rename file', callback: this.renameFile.bind(this), role: 'editor' },
        { label: 'Protect file', callback: this.protectFile.bind(this), role: 'administrator' },
        { label: 'Remove file', callback: this.removeFile.bind(this), role: 'editor' },
      ]
    );
    return <MenuControl items={items}/>;
  }
}

export default connect(
  null,
  (dispatch, {bucket_name}) => ({
    removeFile: (file) => dispatch({ type: 'REMOVE_FILE', file, bucket_name }),
    unprotectFile: (file) => dispatch({ type: 'UNPROTECT_FILE', file, bucket_name }),
    protectFile: (file) => dispatch({ type: 'PROTECT_FILE', file, bucket_name }),
    renameFile: (file, new_file_path) => dispatch({ type: 'RENAME_FILE', file, new_file_path, bucket_name })
  })
)(FileControl);

