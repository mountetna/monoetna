import * as React from 'react';
import { connect } from 'react-redux';
import { copyText } from '../../utils/copy';
import { filePath } from '../../utils/file';
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
        { label: 'Unprotect folder', callback: this.unprotectFolder.bind(this) }
      ] : [
        { label: 'Rename folder', callback: this.renameFolder.bind(this) },
        { label: 'Protect folder', callback: this.protectFolder.bind(this) },
        { label: 'Remove folder', callback: this.removeFolder.bind(this) },
      ];
    return <MenuControl items={items}/>;
  }
}

export default connect(
  null,
  (dispatch) => ({
    removeFolder: (folder) => dispatch({ type: 'REMOVE_FOLDER', folder }),
    unprotectFolder: (folder) => dispatch({ type: 'UNPROTECT_FOLDER', folder }),
    protectFolder: (folder) => dispatch({ type: 'PROTECT_FOLDER', folder }),
    renameFolder: (folder, new_folder_path) => dispatch({ type: 'RENAME_FOLDER', folder, new_folder_path })
  })
)(FolderControl);

