import React, {useContext, useEffect, useState} from 'react';

import markdown from 'etna-js/utils/markdown';

import {VulcanContext} from '../../contexts/vulcan_context';
import {workflowByName} from '../../selectors/workflow_selectors';

export default function Vignette({workflowName}: {workflowName: string}) {
  let {state} = useContext(VulcanContext);
  const workflow = workflowByName(workflowName, state);
  const [text, setText] = useState('');

  useEffect(() => {
    if (workflow) {
      setText(
        workflow.vignette
          ? workflow.vignette
          : 'No vignette provided.'
      );
    }
  }, [workflow, setText]);

  return (
    <div
      className='markdown'
      dangerouslySetInnerHTML={{__html: markdown(text)}}
      style={{
        height: '400px',
        overflowY: 'auto'
      }}
    />
  );
}
