import * as React from 'react';
import { connect } from 'react-redux';
import { selectCurrentFolder } from '../selectors/directory-selector';
import { FolderLink, RootFolderLink } from './folder-link';

const folderCrumb = (folders, folder_name, i) => {
  let folder_path = folders.slice(0,i).join('/');
  let not_last = i < folders.length - 1;

  return <div key={i} className='folder-link' title={folder_name}>
    {
      not_last ?
      <FolderLink folder_name={folder_name} folder_path={folder_path}/>
      :
      <span>{folder_name}</span>
    }
    {
      not_last && <span className='folder-sep fas fa-chevron-right'/>
    }
  </div>;
}

class FolderBreadcrumb extends React.Component {
  render() {
    let { current_folder } = this.props;

    let folders = (current_folder || '').split(/\//).filter(_=>_);
    let root_name = CONFIG.project_name;

    return (
      <div id='folder-breadcrumb'>
        <div key='root' className='folder-link' title={root_name}>
          <div className='folder-project'>
            project
            <span className='folder-sep fas fa-chevron-right'/>
          </div>
          {
            folders.length > 0 ?
            <RootFolderLink name={ root_name }/>
            :
            <span>{root_name}</span>
          }
          {
            folders.length > 0 && <span className='folder-sep fas fa-chevron-right'/>
          }
        </div>
        {
          folders.map((folder_name,i)=> folderCrumb(folders, folder_name, i))
        }
      </div>
    );
  }
}

export default connect(
  // map state
  (state) => ({current_folder: selectCurrentFolder(state)})

)(FolderBreadcrumb);
