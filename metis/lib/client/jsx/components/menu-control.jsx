import * as React from 'react';
import { connect } from 'react-redux';
import { selectUserRole } from '../selectors/user-selector';

const MenuList  = ({items, control, dismiss}) => {
  let { top, right } = control.current.getBoundingClientRect();
  let { documentElement } = document;
  let className = `menu-list ${
    (top + 200 < documentElement.clientHeight) ? 'up' : 'down'
  }`;

  return <ul className={ className } >
    {
      items.map((item,i)=>
        <li onClick={ () => { item.callback(); dismiss(); } } key={i}>
          {item.label}
        </li>
      )
    }
  </ul>;
};

class MenuControl extends React.Component{
  constructor(props) {
    super(props);
    this.control = React.createRef();
    this.state = {
      list_open: false
    };

    this.close = this.close.bind(this);
  }

  componentDidUpdate() {
    let { list_open } = this.state;
    setTimeout(() => {
      if (list_open) window.addEventListener('click', this.close);
      else window.removeEventListener('click', this.close);
    }, 0);
  }

  componentWillUnmount() {
    window.removeEventListener('click', this.close);
  }

  close() { this.setState({ list_open: false }); }

  showDialog() {
    let { items } = this.props;

    this.setState({ list_open: true });
  }

  render() {
    let { items, user_role } = this.props;
    let { list_open } = this.state;
    let className = list_open ? 'control selected' : 'control';

    items = items.filter(({role: item_role}) => !item_role || user_role <= item_role);

    if (!items.length) return null;

    return (
      <div className='control-btn-group'>
        <div className={ className } ref={this.control}
          onClick={this.showDialog.bind(this)}>
          &bull;
          &bull;
          &bull;
        </div>
        {
          list_open && <MenuList control={ this.control } items={items} dismiss={ this.close } />
        }
      </div>
    );
  }
}

export default connect(
  state => ({ user_role: selectUserRole(state) })
)(MenuControl);
