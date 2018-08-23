import * as React from 'react';
import { connect } from 'react-redux';

class ListDialog extends React.Component {
  render() {
    let { items } = this.props;

    return <ul className='list-dialog'>
      {
        items.map((item,i)=> <li onClick={item.callback} key={i}>
          {item.label}
        </li>)
      }
    </ul>
  }
}

class ModalDialog extends React.Component {
  render() {
    let { dismissDialog, dialog, width=150 } = this.props;

    if (!dialog || !Object.keys(dialog).length) return null;

    let { left, top, type, ...dialog_props } = dialog;

    let DialogComponent;

    switch(type) {
      case 'list':
        DialogComponent = ListDialog;
        break;
    }

    let bounds = { top, width };
    let { documentElement } = document;
    console.log(`Document width is ${documentElement.clientWidth}`);
    console.log(documentElement);
    console.log(left);
    console.log(width);
    if (left + width > documentElement.clientWidth) left = left - width;

    if (!DialogComponent) return null;

    return <div className='modal-window' onClick={ dismissDialog }>
      <div className='modal-dialog' onClick={ (e) => e.stopPropagation() }
        style={ { left, top, width: `${width}px` } }>
        <DialogComponent { ...dialog_props }/>
      </div>
    </div>;
  }
}

export default connect(
  ({dialog}) => ({dialog}),
  (dispatch) => ({
    dismissDialog: () => dispatch({type:'DISMISS_DIALOG'})
  })
)(ModalDialog);
