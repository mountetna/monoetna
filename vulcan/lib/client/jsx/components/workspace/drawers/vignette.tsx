import React, {useContext, useEffect, useState} from 'react';

import markdown from 'etna-js/utils/markdown';

import {VulcanContext} from '../../../contexts/vulcan_context';

export default function Vignette({}) {
  let {state} = useContext(VulcanContext);
  const [text, setText] = useState('');

  useEffect(() => {
    if ('vignette.md' in state.status.file_contents) {
      setText(
        state.status.file_contents['vignette.md']
      );
    }
  }, [state.status.file_contents]);

  return (
    <div
      className='markdown'
      dangerouslySetInnerHTML={{__html: markdown(text)}}
      style={{
        maxHeight: '65vh',
        overflowY: 'auto'
      }}
    />
  );
}
