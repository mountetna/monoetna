import React from 'react';
import markdown from 'etna-js/utils/markdown';

export default function MarkdownOutput({data}) {
  return <React.Fragment>
    {Object.keys(data).map(k => {
      const d = data[k];
      if (!d) return null;
      return <div
        key={k}
        // className='markdown'
        style={{paddingLeft: '20px'}}
        dangerouslySetInnerHTML={{__html: markdown(d)}}
        />
    })}
  </React.Fragment>;
}

export function CollapsibleMarkdownOutput({data}) {
  const redata = Object.fromEntries(Object.keys(data).map(k => {
    const d = data[k];
    if (!d) return [k, null];
    return [k, `<details open>\n\n<summary>${k}</summary>\n\n${markdown(d)}\n\n</details>`];
  }))

  return <MarkdownOutput data={redata}/>
}
