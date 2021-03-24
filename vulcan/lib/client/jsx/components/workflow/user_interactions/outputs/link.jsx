import React from 'react';
import Link from 'etna-js/components/link';

export default function LinkOutput({data}) {
  return <Link link={data}>Download data here</Link>;
}
