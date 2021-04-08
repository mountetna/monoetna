import * as React from 'react';
import { connect } from 'react-redux';
import { copyText } from 'etna-js/utils/copy';
import { filePath } from 'etna-js/utils/file';
import MenuControl from '../menu-control';

class FolderControl extends React.Component{
  unprotectFolder() {
    let { unprotectFolder, folder} = this.props;

    unprotectFolder(folder);
  }

  protectFolder() {
    let { protectFolder, folder} = this.props;

    protectFolder(folder);
  }

  renameFolder() {
    let { renameFolder, current_folder, folder } = this.props;
    let new_folder_name = prompt("What is the new name of this folder?", folder.folder_name);
    if (new_folder_name) renameFolder(folder, filePath(current_folder, new_folder_name));
  }

  removeFolder() {
    let { removeFolder, folder} = this.props;

    removeFolder(folder);
  }

  copyLink() {
    let { folder: { download_url } } = this.props;

    copyText(download_url);
  }

  render() {
    let { folder } = this.props;
    let items = folder.read_only ?
      [
        { label: 'Unprotect folder', callback: this.unprotectFolder.bind(this), role: 'administrator' }
      ] : [
        { label: 'Rename folder', callback: this.renameFolder.bind(this), role: 'editor' },
        { label: 'Protect folder', callback: this.protectFolder.bind(this), role: 'administrator' },
        { label: 'Remove folder', callback: this.removeFolder.bind(this), role: 'editor' },
      ];
    return <MenuControl items={items}/>;
  }
}

export default connect(
  null,
  (dispatch, {bucket_name}) => ({
    removeFolder: (folder) => dispatch({ type: 'REMOVE_FOLDER', folder, bucket_name }),
    unprotectFolder: (folder) => dispatch({ type: 'UNPROTECT_FOLDER', folder, bucket_name }),
    protectFolder: (folder) => dispatch({ type: 'PROTECT_FOLDER', folder, bucket_name }),
    renameFolder: (folder, new_folder_path) => dispatch({ type: 'RENAME_FOLDER', folder, new_folder_path, bucket_name })
  })
)(FolderControl);

