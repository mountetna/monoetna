import * as React from 'react';
import * as ReactRedux from 'react-redux';

const ListColumnHead = ({ columnName, width, showName=true }) => {
  let columnLabel = columnName.replace(/\b\w/g, l => l.toUpperCase());
  let columnId = `list-${columnName}-column`;

  return <div id={ columnId } className='list-head-title'
    style={ { flexBasis: width } } >
    { showName && columnLabel }
    { showName && <div className='list-column-head-arrow-group'>
      <span className='fa fa-chevron-down'></span>
    </div>
    }
  </div>
};

export default class ListHead extends React.Component{
  render() {
    let { columns } = this.props;
    return (
      <div id='list-head-group'>
        {
          columns.map( ({name, width, hide}) =>
            <ListColumnHead key={name} width={ width } columnName={ name } showName={ !hide } />
          )
        }
      </div>
    );
  }
}
