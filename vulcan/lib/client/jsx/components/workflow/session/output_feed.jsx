import React, {useContext} from 'react';

import {VulcanContext} from '../../../contexts/vulcan';
import {validPath} from '../../../utils/workflow';
import {completedUiOutputsSelector} from '../../../selectors/workflow';

import StepOutput from '../steps/step_output';
import OutputFocus from './output_focus';

export default function OutputFeed() {
  // Shows stream of Output, Plots, etc.,
  //   as the session object updates.
  const context = useContext(VulcanContext);
  const {workflow, session, pathIndex, status} = context;

  if (!workflow || !validPath({workflow, pathIndex}) || !session || !status)
    return null;

  let outputs = completedUiOutputsSelector(context);

  return (
    <div className='session-output-feed'>
      {outputs.map((s, index) => (
        <StepOutput key={index} step={s.step} stepIndex={s.index}></StepOutput>
      ))}
      <OutputFocus />
    </div>
  );
}
