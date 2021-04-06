import React from 'react';

function stringify(d) {
  try {
    return JSON.stringify(d);
  } catch {
    return d + "";
  }
}

export default function RawOutput({data}) {
  return <React.Fragment>
    {Object.keys(data).map(k => {
      const d = data[k];
      if (!d) return null;

      return <div className='raw-view' key={k}>{stringify(d)}</div>;
    })}
  </React.Fragment>
}
