import React from 'react';
import ConsignmentTable from 'etna-js/plots/components/consignment/consignment_table';
import Consignment from 'etna-js/plots/models/consignment';

export default function ConsignmentOutput({data}) {
  return <React.Fragment>
    {Object.keys(data).map(k => {
      const d = data[k];
      if (!d) return null;

      return (
        <div className='consignment-view' key={k}>
          <ConsignmentTable consignment={new Consignment(d)}/>
        </div>
      );
    })}
  </React.Fragment>
}
