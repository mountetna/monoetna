import * as React from 'react';
import ConsignmentResult from './consignment_result';

export default function ConsignmentView({consignment}) {
  if (!consignment) return null;

  return (
    <React.Fragment>
      {Object.keys(consignment).map((name, i) => (
        <div key={i} className='consignment-variable-group'>
          <div className='consignment-variable-name'>{name}</div>
          <div className='consignment-variable-result'>
            <ConsignmentResult name={name} data={consignment[name]} />
          </div>
        </div>
      ))}
    </React.Fragment>
  );
}
