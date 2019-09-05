import * as React from 'react';
import { connect } from 'react-redux';
import ConfigureBucket from './dialogs/configure-bucket-dialog';
import { camelCase, capitalize } from '../utils/format';

const DIALOGS = { ConfigureBucket };

const ModalDialog = ({ dialog, dismissDialog }) => {
  if (!dialog || !Object.keys(dialog).length) return null;

  let { type, ...dialog_props } = dialog;

  let DialogComponent = DIALOGS[ capitalize(camelCase(type)) ];

  if (!DialogComponent) return null;

  return <div className='modal-window' onClick={ () => dismissDialog() }>
    <div className='modal-dialog' onClick={ (e) => e.stopPropagation() }>
      <DialogComponent { ...dialog_props } />
    </div>
  </div>;
}

export default connect(
  ({dialog}) => ({dialog}),
  (dispatch) => ({
    dismissDialog: () => dispatch({type:'DISMISS_DIALOG'})
  })
)(ModalDialog);
