import * as React from 'react';
import { connect } from 'react-redux';

import ListUpload from './list-upload';
import ListFolder from './list-folder';
import ListFile from './list-file';
import ListUploadFailed from './list-upload-failed';

import { selectFiles, selectFolders, selectUploads } from '../../selectors/directory-selector';

class ListBody extends React.Component{
  render() {
    let { widths, uploads, files, folders, bucket_name, folder_name } = this.props;

    let order = (key) => (a,b) => a[key].localeCompare(b[key]);
    let download_files = Object.values(files).sort(order('file_name'));
    let download_folders = Object.values(folders).sort(order('folder_name'));

    return (
      <div id='list-body-group'>
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
                current_folder={folder_name}
                bucket_name={bucket_name}
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
                current_folder={folder_name}
                bucket_name={bucket_name}
                widths={widths} />
            )
            : null
        }
        { !download_files.length && !download_folders.length && !Object.values(uploads).length && <div className='empty'>No contents</div> }
      </div>
    );
  }
}


export default connect(
  // map state
  (state) => ({
    files: selectFiles(state),
    uploads: selectUploads(state),
    folders: selectFolders(state)
  })
)(ListBody);
