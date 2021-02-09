import React from 'react';

import Plot from './output_renderers/plot';
import Raw from './output_renderers/raw';

export default function Output({viewOption}) {
  if (!viewOption) return null;

  const viewOptions = {
    plot: Plot,
    raw: Raw
  };

  let Component = viewOptions[viewOption];

  return <Component md5sum='some-cell-hash'></Component>;
}
