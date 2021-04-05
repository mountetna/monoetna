import React from 'react';
import Link from 'etna-js/components/link';

export default function LinkOutput({data}) {
  return <React.Fragment>
    {Object.values(data).map(url =>
        <React.Fragment key={url}>
          <Link link={url}>{url.split('/').slice(-1)}</Link>
        </React.Fragment>
    )}
  </React.Fragment>
}
