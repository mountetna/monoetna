import React, {useContext} from 'react';

import {VulcanContext} from '../../../contexts/vulcan_context';
import {dontDownloadForOutputTypes, OUTPUT_TYPES, OUTPUTS} from '../../ui_components';

import StepName from './step_name';

import {statusOfStep, uiComponentOfStep, stepInputDataUrls, stepInputDataRaw} from '../../../selectors/workflow_selectors';
import {WorkspaceStep} from '../../../api_types';

export default function StepOutput({step}: {step: WorkspaceStep}) {
  let {state} = useContext(VulcanContext);
  const vulcan_config = state.workspace?.vulcan_config || {};
  const stepStatus = statusOfStep(step, state.status);
  const uiOutput = uiComponentOfStep(step.name, vulcan_config);

  if (!stepStatus || !uiOutput) return null;
  const stepType = uiOutput in OUTPUTS ? uiOutput : 'default';

  let data;
  if (dontDownloadForOutputTypes.includes(stepType)) {
    data = stepInputDataUrls(step, state.status);
  } else {
    data = stepInputDataRaw(step, state.status, state.data, state.session);
  }

  let url;
  if ([OUTPUT_TYPES.PLOT, OUTPUT_TYPES.PNG].includes(stepType)) {
    url = stepInputDataUrls(step, state.status);
  }

  let OutputComponent = OUTPUTS[stepType];

  return (
    <div className='step-output'>
      <StepName step={step} />
      <div className='outputs-pane'>
        <OutputComponent data={data} url={url}/>
      </div>
    </div>
  );
}
