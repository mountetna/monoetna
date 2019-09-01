import * as React from 'react';
import { connect } from 'react-redux';

class MenuControl extends React.Component{
  constructor(props) {
    super(props);
    this.control = React.createRef();
  }

  showDialog() {
    console.log(this.control);
    let { showDialog, items } = this.props;
    let { top, right } = this.control.current.getBoundingClientRect();

    let dialog = {
      type: 'list',
      items,
      width: 175,
      top, left: right
    };

    showDialog(dialog);
  }

  render() {
    let className = 'control-btn-group';
    return (
      <div ref={this.control}
        className={className}
        onClick={this.showDialog.bind(this)}>
        &bull;
        &bull;
        &bull;
      </div>
    )
  }
}

export default connect(
  null,
  (dispatch) => ({
    showDialog: (dialog) => dispatch({ type: 'SHOW_DIALOG', dialog}),
  })
)(MenuControl);

