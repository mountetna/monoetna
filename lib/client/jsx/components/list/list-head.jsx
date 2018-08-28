import * as React from 'react';
import * as ReactRedux from 'react-redux';

const ListColumnHead = ({ columnName, widths, showName=true }) => {
  let columnLabel = columnName.replace(/\b\w/g, l => l.toUpperCase());
  let columnId = `list-${columnName}-column`;

  return <div id={ columnId } className='list-head-title'
    style={ { flexBasis: widths[columnName] } } >
    { showName && columnLabel }
    { showName && <div className='list-column-head-arrow-group'>
      <span className='fa fa-chevron-down'></span>
    </div>
    }
  </div>
};

export default class ListHead extends React.Component{
  render() {
    let { widths } = this.props;
    return (
      <div id='list-head-group'>
        <ListColumnHead widths={ widths } columnName='type'/>
        <ListColumnHead widths={ widths } columnName='name'/>
        <ListColumnHead widths={ widths } columnName='status' showName={false}/>
        <ListColumnHead widths={ widths } columnName='updated'/>
        <ListColumnHead widths={ widths } columnName='size'/>
        <ListColumnHead widths={ widths } columnName='control' showName={false}/>
      </div>
    );
  }
}
