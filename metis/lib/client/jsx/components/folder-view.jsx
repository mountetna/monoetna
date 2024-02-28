import React, {useEffect, useCallback, useRef, useState} from 'react';
import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
import {useReduxState} from 'etna-js/hooks/useReduxState';
import {useAsync} from 'etna-js/utils/cancellable_helpers'

import ListBody from './list/list-body';
import ListHead from './list/list-head';
import FolderBreadcrumb from './folder-breadcrumb';
import ControlBar from './control-bar';

import {selectCurrentFolder} from '../selectors/directory-selector';
import AutorenewIcon from '@material-ui/icons/Autorenew';
import { makeStyles } from '@material-ui/core/styles';

const useStyles = makeStyles((theme) => ({
  loadingIcon: {
    'animation': '$spin 4s linear infinite'
  },
  '@keyframes spin': {
      '100%': {
          '-webkit-transform': 'rotate(360deg)',
          'transform': 'rotate(360deg)',
      }
  },
}))

const COLUMNS = [
  {name: 'type', width: '60px'},
  {name: 'name', width: '60%'},
  {name: 'status', width: '90px', hide: true},
  {name: 'updated', width: '30%'},
  {name: 'size', width: '10%'},
  {name: 'control', width: '100px', hide: true}
];

const COLUMN_WIDTHS = COLUMNS.reduce((widths, column) => {
  widths[column.name] = column.width;
  return widths;
}, {});

const INVALID = '\ninvalid\n';

const InvalidFolder = () => (
  <div className='invalid-folder-view-group'>Invalid folder!</div>
);

const FolderView = ({bucket_name, folder_name}) => {
  const invoke = useActionInvoker();
  const current_folder = useReduxState((state) => selectCurrentFolder(state));
  const uploadFileInput = useRef(null);
  const uploadDirInput = useRef(null);

  const [loading, setLoading] = useState(true)
  const [showLoading, setShowLoading] = useState(false)

  const classes = useStyles();

  useEffect(() => {
    setLoading(true);
    invoke({type: 'RETRIEVE_FILES', bucket_name, folder_name});
    setLoading(false);
  }, []);

  useAsync( async () => {
    if (loading && !showLoading) {
      const delay = (delayInms) => {
        return new Promise(resolve => setTimeout(resolve, delayInms));
      };
      let delayres = await delay(100);
      setShowLoading(true);
    }
    if (!loading) {
      setShowLoading(false);
    }
  }, [loading, showLoading]);

  const selectUpload = useCallback(() => {
    invoke({
      type: 'SHOW_DIALOG',
      dialog: {
        type: 'upload_dialog',
        startDirectoryUpload: () => uploadDirInput.current.click(),
        startFileUpload: () => uploadFileInput.current.click()
      }
    });
  }, [invoke]);

  const prepareFiles = useCallback(
    (event, input) => {
      if (event === undefined) return;
      let {files} = input.current;

      for (let i = 0; i < files.length; i++) {
        invoke({
          type: 'FILE_SELECTED',
          file: files[i],
          folder_name,
          bucket_name
        });
      }

      // Reset the input field.
      input.value = '';
    },
    [invoke, folder_name, bucket_name]
  );

  const fileSelected = useCallback(
    (event) => {
      prepareFiles(event, uploadFileInput);
    },
    [prepareFiles]
  );

  const dirSelected = useCallback(
    (event) => {
      prepareFiles(event, uploadDirInput);
    },
    [prepareFiles]
  );

  const selectFolder = useCallback(() => {
    let new_folder_name = prompt('Enter the folder name', 'Untitled Folder');

    if (new_folder_name) {
      invoke({
        type: 'CREATE_FOLDER',
        parent_folder: folder_name,
        folder_name: new_folder_name,
        bucket_name
      });
    }
  }, [invoke, folder_name, bucket_name]);

  const selectFolderDownload = useCallback(() => {
    invoke({type: 'LIST_FILES_RECURSIVE', folder_name, bucket_name}).then(
      (files) => {
        invoke({type: 'DOWNLOAD_FILES_ZIP', files, folder_name, bucket_name});
      }
    );
  }, [invoke, folder_name, bucket_name]);

  const copyFile = useCallback(() => {
    let metis_path = prompt(
      'Enter the Metis path of the file to paste in this folder'
    );

    if (metis_path) {
      invoke({
        type: 'COPY_FILE',
        folder_name,
        bucket_name,
        metis_path,
        dest_metis_path: `metis://${
          CONFIG.project_name
        }/${bucket_name}/${folder_name}/${metis_path.split('/').at(-1)}`
      });
    }
  }, [invoke, folder_name, bucket_name]);

  if (current_folder == INVALID) return <InvalidFolder />;

  let buttons = [
    {
      onClick: selectFolder,
      title: 'Create folder',
      icon: 'folder',
      overlay: 'plus',
      role: 'editor'
    },
    {
      onClick: copyFile,
      title: 'Paste file with Metis path',
      icon: 'cloud',
      overlay: 'paste',
      role: 'editor'
    },
    {
      onClick: selectUpload,
      title: 'Upload file(s)',
      icon: 'upload',
      role: 'editor'
    },
    {
      onClick: selectFolderDownload,
      title: 'Download directory as zip',
      icon: 'download',
      role: 'viewer'
    }
  ];

  return (
    <div className='folder-view-group'>
      <div className='control-group'>
        <FolderBreadcrumb folder_name={folder_name} bucket_name={bucket_name} />
        <ControlBar buttons={buttons}>
          {/* For uploading individual files */}
          <input
            name='upload-file'
            type='file'
            multiple='multiple'
            style={{display: 'none'}}
            ref={uploadFileInput}
            onChange={fileSelected}
          />

          {/* For uploading directories */}
          <input
            name='upload-directory'
            type='file'
            webkitdirectory='webkitdirectory'
            directory='directory'
            multiple='multiple'
            style={{display: 'none'}}
            ref={uploadDirInput}
            onChange={dirSelected}
          />
        </ControlBar>
      </div>
      <div className='listing-group'>
        <ListHead columns={COLUMNS} />
        <ListBody
          widths={COLUMN_WIDTHS}
          folder_name={folder_name}
          bucket_name={bucket_name}
        />
        {loading ? <>
          <AutorenewIcon className={classes.loadingIcon}/>
          Loading
        </> : null}
      </div>
    </div>
  );
};

export default FolderView;
