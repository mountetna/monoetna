import * as React from 'react';
import { connect } from 'react-redux';
import { FolderLink, RootFolderLink } from './folder-link';

const Chevron = () => <span className='folder-sep fas fa-chevron-right'/>;

const FolderCrumb = ({folder_path, link, bucket_name, folder_name}) =>
  <div className='folder-link' title={folder_name}>
    {
      link
      ? <FolderLink log={true} bucket_name={bucket_name} folder_name={folder_name} folder_path={folder_path}/>
      : <span>{folder_name || bucket_name}</span>
    }
    {
      link && <Chevron/>
    }
  </div>;

const FolderBreadcrumb = ({ folder_name, bucket_name }) => {
  let folders = (folder_name || '').split(/\//).filter(_=>_);
  let root_name = CONFIG.project_name;

  return (
    <div id='folder-breadcrumb'>
      <div key='root' className='folder-link' title={root_name}>
        <div className='folder-project'>
          project
          <Chevron/>
        </div>
        {
          bucket_name
            ? <RootFolderLink name={ root_name }/>
            : <span>{ root_name }</span>
        }
        { bucket_name && <Chevron/> }
      </div>
      <FolderCrumb link={ folders.length > 0 } bucket_name={ bucket_name } />
      {
        folders.map((folder_name,i)=> <FolderCrumb
          key={i}
          folder_path={ folders.slice(0,i).join('/') }
          link={ i < folders.length - 1}
          bucket_name={ bucket_name }
          folder_name={ folder_name }
          />)
      }
    </div>
  );
};

export default FolderBreadcrumb;
