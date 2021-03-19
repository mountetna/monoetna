import React from 'react';

export default function RawOutput({data}) {
  if (!data) return null;
  return <div className='raw-view'>{data}</div>;
}
