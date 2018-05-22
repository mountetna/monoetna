import * as React from 'react';
import { connect } from 'react-redux';

const FolderLink = (folders, i) => {
  let folder_path = folders.slice(1,i+1).join('/');
  let folder_name = folders[i];

  let folder_url = `/${CONFIG.project_name}/${i ? 'browse/' : ''}${folder_path}`;

  return <div key={i} className='folder-link' title={folder_name}>
    {
      (i < folders.length-1) ?
      <a href={ folder_url }>{folder_name}</a>
      :
      <span>{folder_name }</span>
    }
    {
      (i < folders.length-1) && <span className='folder-sep fas fa-chevron-right'/>
    }
  </div>;
}

class FolderBreadcrumb extends React.Component {
  render() {
    let { current_folder } = this.props;

    let folders = (current_folder || '').split(/\//).filter(_=>_);

    folders.unshift('All files');

    console.log(folders);

    return (
      <div id='folder-breadcrumb'>
        {
          folders.map((folder_name,i)=>
            FolderLink(folders, i)
          )
        }
      </div>
    );
  }
}

export default connect(
  // map state
  ({files: { current_folder }}) => ({current_folder})

)(FolderBreadcrumb);
