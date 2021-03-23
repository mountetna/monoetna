import React from 'react';
import Link from 'etna-js/components/link';

export default function LinkOutput({data}) {
  return <React.Fragment>
    {Object.values(data).map(url =>
        <React.Fragment key={url}>
          <Link link={url}>Download data here</Link>
        </React.Fragment>
    )}
  </React.Fragment>
}
