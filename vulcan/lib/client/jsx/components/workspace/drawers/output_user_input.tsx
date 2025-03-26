import React, {useContext} from 'react';

import {VulcanContext} from '../../../contexts/vulcan_context';
import {dataAndUrlForOutputTypes, dontDownloadForOutputTypes, OUTPUTS} from '../../ui_components';

import StepIconName from './step_elements/step_icon_name';

import {statusOfStep, uiComponentOfStep, stepInputDataUrls, stepInputDataRaw} from '../../../selectors/workflow_selectors';
import {WorkspaceStep} from '../../../api_types';

export default function OutputUI({step}: {step: WorkspaceStep}) {
  let {state} = useContext(VulcanContext);
  let {workspace, status} = state;
  if (!workspace?.vulcan_config) return null;
  const vulcan_config = workspace.vulcan_config;
  const uiOutput = uiComponentOfStep(step.name, vulcan_config);

  if (!uiOutput) return null;
  const stepType = uiOutput in OUTPUTS ? uiOutput : 'default';

  let data;
  if (dontDownloadForOutputTypes.includes(stepType)) {
    data = step.input.files;
  } else {
    data = stepInputDataRaw(step, status.last_params, status.file_contents);
  }

  let url;
  if (dataAndUrlForOutputTypes.includes(stepType)) {
    url = step.input.files;
  }

  let OutputComponent = OUTPUTS[stepType];

  return (
    <div className='step-output'>
      <StepIconName step={step} />
      <div className='outputs-pane'>
        <OutputComponent data={data} url={url}/>
      </div>
    </div>
  );
}
