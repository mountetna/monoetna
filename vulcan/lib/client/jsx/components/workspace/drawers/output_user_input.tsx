import React, {useContext} from 'react';

import {VulcanContext} from '../../../contexts/vulcan_context';
import {dataAndUrlForOutputTypes, dontDownloadForOutputTypes, OUTPUTS} from '../../ui_components';

import StepIconName from './step_elements/step_icon_name';

import {configIOValues, uiComponentOfStep} from '../../../selectors/workflow_selectors';
import {VulcanConfigElement, WorkspaceStep} from '../../../api_types';
import { DataEnvelope, fillInputData } from '../ui_definitions/input_types';

export default function OutputUI({step}: {step: WorkspaceStep}) {
  const {state} = useContext(VulcanContext);
  const {workspace, status} = state;
  const config = workspace?.vulcan_config ? workspace.vulcan_config : undefined;

  if (!config) return null;

  const uiOutput = uiComponentOfStep(step.name, config);

  if (!uiOutput) return null;
  const stepType = uiOutput in OUTPUTS ? uiOutput : 'default';

  let data: DataEnvelope<any>;
  if (dontDownloadForOutputTypes.includes(stepType)) {
    data = configIOValues(step.input.files);
  } else {
    data = fillInputData(step as VulcanConfigElement, step, status.last_params, status.file_contents);
  }

  let url;
  if (dataAndUrlForOutputTypes.includes(stepType)) {
    url = configIOValues(step.input.files);
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
