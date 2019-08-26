import * as React from 'react';
import { connect } from 'react-redux';

class ListDialog extends React.Component {
  clickItem({callback, show}) {
    callback();
    if (!show) this.props.dismissDialog();
  }

  render() {
    let { items, onClick } = this.props;
    return <ul className='list-dialog'>
      {
        items.map((item,i)=> <li onClick={() => this.clickItem(item) } key={i}>
          {item.label}
        </li>)
      }
    </ul>;

  }
}

export default connect(
  null,
  (dispatch) => ({
    dismissDialog: () => dispatch({type:'DISMISS_DIALOG'})
  })
)(ListDialog);
