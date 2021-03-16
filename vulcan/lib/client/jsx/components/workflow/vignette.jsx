import React, {useContext, useEffect, useState} from 'react';

import markdown from 'etna-js/utils/markdown';

import {VulcanContext} from '../../contexts/vulcan';
import {workflowByName} from '../../selectors/workflow';

export default function Vignette({workflowName}) {
  let {workflows} = useContext(VulcanContext);

  const [text, setText] = useState('No vignette provided.');

  useEffect(() => {
    if (
      workflows.workflows &&
      workflowByName({
        workflows: workflows.workflows,
        workflowName
      })
    ) {
      let selectedWorkflow = workflowByName({
        workflows: workflows.workflows,
        workflowName
      });
      setText(
        selectedWorkflow.vignette
          ? selectedWorkflow.vignette
          : 'No vignette provided.'
      );
    }
  }, [workflows]);

  return (
    <div
      className='markdown'
      dangerouslySetInnerHTML={{__html: markdown(text)}}
    />
  );
}
