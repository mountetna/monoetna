import * as React from 'react';
import { connect } from 'react-redux';

class MenuControl extends React.Component{
  constructor(props) {
    super(props);
    this.state = { dialog_shown: false };
  }

  dismissDialog() {
    this.setState({dialog_shown: false});
  }

  showDialog() {
    let { showDialog, items } = this.props;
    let { top, left } = this.control.getBoundingClientRect();

    this.setState({ dialog_shown: true });

    let dialog = {
      type: 'list',
      items,
      dismiss: this.dismissDialog.bind(this),
      top,
      left
    };

    showDialog(dialog);
  }

  setControlRef(ref) {
    this.control = ref;
  }

  render() {
    let { dialog_shown } = this.state;
    let className = `control-btn-group ${dialog_shown ? 'selected' : ''}`;
    return (
      <div ref={ this.setControlRef.bind(this) } className={className} onClick={this.showDialog.bind(this)}>
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

