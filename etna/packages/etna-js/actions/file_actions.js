import React, {useRef, useState, useCallback, useEffect} from 'react';
import Modal from 'react-modal';

import ButtonBar from '../components/button_bar';
import Icon from '../components/icon';

export const STUB = '::blank';
export const TEMP = '::temp';
import {PickBucket, PickFileOrFolder} from '../components/metis_exploration';

// We don't have a lot of content, so let's get a smaller Modal
export const customStyles = {
  content: {
    top: '50%',
    left: '50%',
    right: 'auto',
    bottom: 'auto',
    marginRight: '-50%',
    minWidth: '40%',
    transform: 'translate(-50%, -50%)'
  }
};

// metis:\/\/([^\/]*?)\/([^\/]*?)\/(.*)$
export const METIS_PATH_MATCH = (path) =>
  new RegExp(
    '^metis://' +
      // project_name
      '([^/]*?)/' +
      // bucket_name
      '([^/]*?)/' +
      // folder path + filename
      '(.*)$'
  ).test(path);

export const useFileInputActions = (
  metis,
  error,
  setMetis,
  setError,
  onChange,
  onBlur
) => {
  const [metisPath, setMetisPath] = useState('');

  const [bucketName, setBucketName] = useState('');
  const [path, setPath] = useState('');
  const [targetType, setTargetType] = useState(null);
  
  function changeBucket(e) {
    setPath('');
    setTargetType(null);
    setBucketName(e);
  }
  function reset() {
    changeBucket('');
  }
  
  useEffect(() => {
    if (targetType=='file') {
      setMetisPath(`metis://${CONFIG.project_name}/${bucketName}/${path}`)
    } else {
      setMetisPath('')
    }
  }, [bucketName, path, targetType])

  return {
    metisSelector,
    closeModal,
    selectMetisFile,
    formatFileRevision,
    setTempRevision,
    isTempRevision
  };

  function isTempRevision(revision) {
    if (!(revision && revision.path)) return false;

    return (
      revision.path.indexOf('/upload/') > -1 &&
      revision.path.indexOf('X-Etna-Signature') > -1
    );
  }

  function metisSelector() {
    // TODO: would be nice to make this like a folder / file search
    return (
      <Modal
        isOpen={metis}
        contentLabel='Enter Metis path'
        style={customStyles}
        onRequestClose={closeModal}
        appElement={document.querySelector('#root')}
      >
        <div className='attribute modal file-metis-select'>
          <h2>Select a Metis file</h2>
          <div className='input-box-wrapper'>
            <PickBucket
              bucket={bucketName}
              label="Bucket"
              setBucket={(e) => changeBucket(e)}
            />
            <PickFileOrFolder
              bucket={bucketName}
              label="Destination Folder"
              useTargetType={setTargetType}
              setPath={(e) => setPath(e)}
              basePath={''}
              topLevelPlaceholder={'top-level of bucket'}
              path={path}
            />
            <div className='modal-button-wrapper'>
              <ButtonBar
                className='modal-buttons'
                buttons={[
                  {type: 'check', click: () => selectMetisFile()},
                  {type: 'cancel', click: () => closeModal()}
                ]}
              />
              {error ? (
                <p className='file-metis-error'>Invalid Metis path</p>
              ) : (
                ''
              )}
              <div className='modal-buttons pull-right'>
                <Icon
                  className=''
                  icon='question-circle'
                  title='Help'
                  onClick={() =>
                    window.open(
                      'https://mountetna.github.io/timur.html#managing-data-files',
                      '_blank'
                    )
                  }
                />
              </div>
            </div>
          </div>
        </div>
      </Modal>
    );
  }

  function closeModal() {
    setMetis(false);
    reset();
  }

  function selectMetisFile() {
    if (!METIS_PATH_MATCH(metisPath)) {
      setError(true);
      return;
    } else {
      setError(false);
      setMetis(false);
    }

    onChange(formatFileRevision(metisPath));
    onBlur();

    reset();
  }

  function formatFileRevision(newValue, files) {
    if (null === newValue) return {path: newValue};

    let file_parts = newValue.split('/');
    let revision = {path: newValue};

    if (STUB !== newValue && TEMP !== newValue && !files) {
      revision['original_filename'] = file_parts[file_parts.length - 1];
    }
    if (files && TEMP === newValue) {
      revision['original_files'] = files;
    }

    return revision;
  }

  function setTempRevision(e) {
    e.preventDefault();

    onChange(formatFileRevision(TEMP, e.target.files));
    onBlur();
  }
};
