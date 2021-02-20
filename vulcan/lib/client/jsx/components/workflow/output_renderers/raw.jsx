import React from 'react';

import ConsignmentView from 'etna-js/plots/components/consignment/consignment_view';

export default function Raw({md5sum}) {
  if (!md5sum) return null;
  return <ConsignmentView md5sum='some-cell-hash'></ConsignmentView>;
}
