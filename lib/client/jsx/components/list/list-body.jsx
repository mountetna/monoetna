import * as React from 'react';
import * as ReactRedux from 'react-redux';

import { ListUpload, ListEntry, ListFolder } from './list-entry';
import ListUploadFailed from './list-upload-failed';

class ListBody extends React.Component{
  constructor() {
    super();
  }

  render() {
    let { files, widths, user } = this.props;
    let { uploads, downloads, fails } = files;
    let { permissions } = user;

    let order = (a,b) => a.file_name.localeCompare(b.file_name);
    let download_files = Object.values(downloads).filter(f=>!f.is_folder).sort(order);
    let download_folders = Object.values(downloads).filter(f=>f.is_folder).sort(order);

    return (
      <div id='list-body-group'>
        {/* Render the incomplete uploads. */}
        { (Object.values(uploads).length) ?
            Object.values(uploads).map( upload =>
              <ListUpload
                key={upload.file_name}
                upload={upload}
                widths={widths}
                user={ user }
                permissions={ permissions } />
            )
            : null
        }
        {
          (download_folders.length) ?
            download_folders.map( folder =>
              <ListFolder
                key={folder.file_name}
                file={folder}
                widths={widths} />
            )
            : null
        }
        {
          (download_files.length) ?
            download_files.map( file =>
              <ListEntry
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


const ListBodyContainer = ReactRedux.connect(
  // map state
  ({user,files}) => ({user,files})
)(ListBody);

export default ListBodyContainer;
