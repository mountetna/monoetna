import React, {useContext} from 'react';

import {VulcanContext} from '../../../contexts/vulcan_context';
import {dataAndUrlForOutputTypes, DataEnvelope, dontDownloadForOutputTypes, OUTPUTS} from '../../ui_components';

import StepIconName from './step_elements/step_icon_name';

import {configIOValues, uiComponentOfStep} from '../../../selectors/workflow_selectors';
import {WorkspaceStatus, WorkspaceStep} from '../../../api_types';

function stepInputDataRaw(
  step: WorkspaceStep,
  last_params: WorkspaceStatus['last_params'],
  file_contents: WorkspaceStatus['file_contents']
): {[k: string]: any} {
  // Pull out any previous step's output data link that is a required
  //   input into this UI step.
  const result: {[k: string]: any} = {};

  configIOValues(step.input.params).forEach((paramName) => {
    result[paramName] = last_params[paramName] || null;
  });

  configIOValues(step.input.files).forEach((fileName) => {
    result[fileName] = file_contents[fileName] || null;
  });

  return result;
};

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
    // const newData = fillInputData(step as VulcanConfigElement, step, status.last_params, status.file_contents);
    const oldData = stepInputDataRaw(step, status.last_params, status.file_contents);
    // console.log({newData})
    // console.log({oldData})
    data = oldData;
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
