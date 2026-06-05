import React, {useContext, useEffect, useState} from 'react';

import markdown from 'etna-js/utils/markdown';

import {VulcanContext} from '../../../contexts/vulcan_context';

export default function Vignette({}) {
  let {state} = useContext(VulcanContext);
  const [text, setText] = useState('');

  useEffect(() => {
    if (!!state.workspace && !!state.workspace.vignette) {
      setText(
        state.workspace.vignette
      );
    }
  }, [state.workspace]);

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
