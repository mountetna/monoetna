import * as React from 'react';
import { connect } from 'react-redux';
import ListDialog from './dialogs/list-dialog';
import ConfiureBucketDialog from './dialogs/configure-bucket-dialog';

const ModalDialog = ({ dialog, dismissDialog }) => {
  if (!dialog || !Object.keys(dialog).length) return null;

  let { left, right, top, bottom, type, height, width=150, ...dialog_props } = dialog;

  let DialogComponent;

  switch(type) {
    case 'list':
      DialogComponent = ListDialog;
      break;
    case 'create-bucket':
      DialogComponent = ConfiureBucketDialog;
      dialog_props.mode = 'create';
      break;
    case 'configure-bucket':
      DialogComponent = ConfiureBucketDialog;
      dialog_props.mode = 'configure';
      break;
  }

  let { documentElement } = document;
  if (left + width > documentElement.clientWidth) left = left - width;

  if (!DialogComponent) return null;

  return <div className='modal-window' onClick={ () => dismissDialog() }>
    <div className='modal-dialog' onClick={ (e) => e.stopPropagation() }
      style={ { left, right, top, bottom, height: height ? `${height}px` : undefined, width: `${width}px` } }>
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
