import React from 'react';

export default function RawOutput({data}) {
  return <React.Fragment>
    {Object.keys(data).map(k => {
      const d = data[k];
      if (!d) return null;

      return <div className='raw-view' key={k}>{data}</div>;
    })}
  </React.Fragment>
}
