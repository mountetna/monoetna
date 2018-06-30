import * as React from 'react';
import { connect } from 'react-redux';

import { ListUpload, ListFile, ListFolder } from './list-entry';
import ListUploadFailed from './list-upload-failed';

import { selectFiles, selectFolders, selectUploads, selectCurrentFolder } from '../../selectors/directory-selector';

class ListBody extends React.Component{
  render() {
    let { widths, uploads, files, folders, current_folder } = this.props;

    let order = (a,b) => a.file_name.localeCompare(b.file_name);
    let download_files = Object.values(files).sort(order);
    let download_folders = Object.values(folders).sort(order);

    return (
      <div id='list-body-group'>
        {/* Render the incomplete uploads. */}
        { (Object.values(uploads).length) ?
            Object.values(uploads).map( upload =>
              <ListUpload
                key={upload.file_name}
                upload={upload}
                widths={widths}
              />
            )
            : null
        }
        {
          (download_folders.length) ?
            download_folders.map( folder =>
              <ListFolder
                key={folder.folder_name}
                folder={folder}
                current_folder={current_folder}
                widths={widths} />
            )
            : null
        }
        {
          (download_files.length) ?
            download_files.map( file =>
              <ListFile
                key={file.file_name}
                file={file}
                widths={widths} />
            )
            : null
        }
      </div>
    );
  }
}


export default connect(
  // map state
  (state) => ({
    files: selectFiles(state),
    uploads: selectUploads(state),
    folders: selectFolders(state),
    current_folder: selectCurrentFolder(state)
  })
)(ListBody);
