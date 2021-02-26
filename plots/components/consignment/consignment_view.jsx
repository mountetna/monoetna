import * as React from 'react';
import { connect } from 'react-redux';
import ConsignmentResult from './consignment_result';
import { selectConsignment } from '../../selectors/consignment_selector';
import Consignment from './consignment'

class ConsignmentView extends React.Component {
  render() {
    let { consignment } = this.props;

    if (!consignment) return null;

    return <div className='consignment-view'>
      <div className='consignment-view-label'>Results:</div>
      <Consignment consignment={consignment}></Consignment>
    </div>
  }
}

export default connect(
  (state, {md5sum}) => ({
    consignment: md5sum && selectConsignment(state, md5sum)
  })
)(ConsignmentView);
