import * as React from 'react';
import { connect } from 'react-redux';
import { selectConsignment } from '../../selectors/consignment_selector';
import ConsignmentTable from './consignment_table'

class ConsignmentView extends React.Component {
  render() {
    let { consignment } = this.props;

    if (!consignment) return null;

    return <div className='consignment-view'>
      <div className='consignment-view-label'>Results:</div>
      <ConsignmentTable consignment={consignment}></ConsignmentTable>
    </div>
  }
}

export default connect(
  (state, {md5sum}) => ({
    consignment: md5sum && selectConsignment(state, md5sum)
  })
)(ConsignmentView);
