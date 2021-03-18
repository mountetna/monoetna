import React from 'react';
import ConsignmentTable from 'etna-js/plots/components/consignment/consignment_table';
import Consignment from 'etna-js/plots/models/consignment';

export default function ConsignmentOutput({data}) {
  if (!data) return null;

  return (
    <div className='consignment-view'>
      <ConsignmentTable consignment={new Consignment(data)}></ConsignmentTable>
    </div>
  );
}
