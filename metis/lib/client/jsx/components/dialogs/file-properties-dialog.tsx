import React, {useCallback} from 'react';

import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';

import {File} from '../../types/metis-types';
import ConfigRow from './config-row';

const FilePropertiesDialog = ({file}: {file: File}) => {
  const invoke = useActionInvoker();

  const close = useCallback(() => {
    invoke({type: 'DISMISS_DIALOG'});
  }, []);

  return (
    <div className='file-properties-dialog'>
      <div className='title'>File properties</div>
      {Object.entries(file).map(([key, value], index) => {
        return (
          <ConfigRow label={key} key={index}>
            <input type='text' readOnly={true} value={value?.toString()} />
          </ConfigRow>
        );
      })}
      <div className='close'>
        <span className={'button'} onClick={close}>
          Close
        </span>
      </div>
    </div>
  );
};

export default FilePropertiesDialog;
