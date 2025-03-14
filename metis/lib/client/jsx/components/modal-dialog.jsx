import * as React from 'react';
import {connect} from 'react-redux';
import {camelCase, capitalize} from '../utils/format';

import ConfigureBucket from './dialogs/configure-bucket-dialog';
import Message from './dialogs/message-dialog';
import UploadDialog from './dialogs/upload-dialog';
import MoveFolder from './dialogs/move-folder-dialog';
import MoveFile from './dialogs/move-file-dialog';
import FileProperties from './dialogs/file-properties-dialog';
import CopyFolder from './dialogs/copy-folder-dialog';

const DIALOGS = {
  ConfigureBucket,
  Message,
  UploadDialog,
  MoveFolder,
  MoveFile,
  FileProperties,
  CopyFolder
};

const ModalDialog = ({dialog, dismissDialog}) => {
  if (!dialog || !Object.keys(dialog).length) return null;

  let {type, ...dialog_props} = dialog;

  let DialogComponent = DIALOGS[capitalize(camelCase(type))];

  if (!DialogComponent) {
    console.error(
      'Modal',
      type,
      'could not be resolved.  Is it in DIALOGS?  Also it needs to be in snake casing'
    );
    return null;
  }

  return (
    <div className='modal-window' onClick={() => dismissDialog()}>
      <div className='modal-dialog' onClick={(e) => e.stopPropagation()}>
        <DialogComponent {...dialog_props} />
      </div>
    </div>
  );
};

export default connect(
  ({dialog}) => ({dialog}),
  (dispatch) => ({
    dismissDialog: () => dispatch({type: 'DISMISS_DIALOG'})
  })
)(ModalDialog);
